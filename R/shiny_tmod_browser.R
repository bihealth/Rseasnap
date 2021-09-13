## call .tmod_browser_prepare_res_single for every data set
.tmod_browser_prepare_res <- function(but, tmod_res) {
  tmod_res <- imap(tmod_res, ~ {
    .tmod_browser_prepare_res_single(.y, but, .x)
  })

  tmod_res
}

## Construct the results table to display. Specifically, add action button
## for launching the plot.
.tmod_browser_prepare_res_single <- function(ds_id, but, tmod_res) {
  # prepare the tmod res
  tmod_res <- tmod_res %>% imap(~ {
    .cntr <- .y
    imap(.x, ~ {
           .dbname <- .y
           imap(.x, ~ {
                  .sort <- .y
                  .x %>% mutate(">"=sprintf(but, ds_id, .data[["ID"]], .cntr, .dbname, .sort)) %>%
                    select(-cES, -cerno) %>%
                    arrange(P.Value) %>%
                    relocate(all_of(">"), .before=1)
            })
         })
  })

  return(tmod_res)
}



#' @rdname tmodBrowserTableServer
#' @export
tmodBrowserTableUI <- function(id, cntr_titles, upset_pane=FALSE) {

  cntr_titles <- .prep_cntr_titles(cntr_titles)

  but <- actionButton("uselessID", label=" \U25B6 ", class = "btn-primary btn-sm")

  if(upset_pane == TRUE) {
    main_pane <-  tabsetPanel(id=NS(id, "main_tabset"),
                       tabPanel("Results", 
                                column(dataTableOutput(NS(id, "tmodResTab")), width=12)),
                       tabPanel("Upset plot", 
                                fluidRow(plotOutput(NS(id, "upset_plot"), height="100%")))
                     )
  } else {
    main_pane <- column(dataTableOutput(NS(id, "tmodResTab")), width=12)
  }


  ui <- sidebarLayout(
          sidebarPanel(
           fluidRow(selectInput(NS(id, "contrast"), label="Contrast", choices=cntr_titles, width="100%")),
           fluidRow(
                    column(uiOutput(NS(id, "table_sel_db")), width=6),
                    column(uiOutput(NS(id, "table_sel_sort")), width=6)
                    ),
           fluidRow(
             checkboxInput(NS(id, "filter"), label="Filter results", value=TRUE),
             fluidRow(
                      column(numericInput(NS(id, "f_auc"),  label="Filter by AUC", 
                          min=.5, max=1.0, step=0.1, value=0.65, width="50%"), width=6),
                      column(numericInput(NS(id, "f_pval"), label="Filter by FDR", 
                          min=0, max=1.0, step=0.1, value=0.05, width="50%"), width=6)
             )
           ),
           HTML(paste("Click on the", as.character(but), "buttons to view an evidence plot")),
           width=3
          ),
          mainPanel(
            fluidRow(main_pane),
            width=9
          )
        )

  return(ui)
}

#' Shiny Module – tmod results browser table selection
#'
#' Shiny Module – tmod results browser table selection
#'
#' @param gs_id a list of reactive values (returned by `reactiveValues()`), including 
#' dataset (`ds`), gene set ID (`id`), contrast id (`cntr`), database ID
#' (`db`) and sorting mode (`sort`). If `mod_id` is not `NULL`, these
#' reactive values will be populated, possibly triggering an action in
#' another shiny module.
#' @param tmod_res results of tmod analysis, returned by `get_tmod_res`
#' @param cntr_titles possibly named character vector with contrast names
#' @param id identifier for the namespace of the module
#' @param tmod_dbs list of lists of tmod database objects (or list of lists
#' of list of tmod database object in multilevel mode). If NULL, upset
#' plots cannot be generated.
#' @param multilevel if TRUE, the results are grouped in data sets
#' @export
tmodBrowserTableServer <- function(id, tmod_res, gs_id=NULL, multilevel=FALSE, tmod_dbs=NULL) {

  if(!multilevel) {
    tmod_res <- list(default=tmod_res)
    tmod_dbs <- list(default=tmod_dbs)
  }

  but <- actionButton("go_%s-!-%s-!-%s-!-%s-!-%s", label=" \U25B6 ", 
                      onclick=sprintf('Shiny.onInputChange(\"%s-select_button\",  this.id)', id),  
                      class = "btn-primary btn-sm")

  tmod_res <- .tmod_browser_prepare_res(as.character(but), tmod_res)

  moduleServer(id, function(input, output, session) {
    message("Launching tmod browser server")

    observeEvent(input$filter, {
                   if(input$filter) {
                     enable("f_auc")
                     enable("f_pval")
                   } else {
                     disable("f_auc")
                     disable("f_pval")
                   }
    })

    dataset  <- reactiveVal()
    contrast <- reactiveVal()

    observeEvent(input$contrast, {
      dataset(gsub("::.*", "", input$contrast))
      contrast(gsub(".*::", "", input$contrast))
    })

    res <- reactiveVal()

    observe({
      if(!(
           isTruthy(dataset()) &&
           isTruthy(contrast()) &&
           isTruthy(input$db) &&
           isTruthy(input$sort)
         )) { return(NULL) }
         
      .res <- tmod_res[[dataset()]][[contrast()]][[input$db]][[input$sort]] 
      
      if(input$filter) {
        .res <- .res %>%
        filter(.data[["AUC"]] > input$f_auc & .data[["adj.P.Val"]] < input$f_pval)
      }
      mf("res are changing, %d gene sets", nrow(.res))

      res(list(
               db=input$db,
               contrast=contrast(),
               sort=input$sort,
               ds=dataset(),
               res=.res))

    })

    output$upset_plot <- renderPlot({
      if(is.null(.res <- res())) { return(NULL) }
      if(is.null(tmod_dbs[[.res$ds]])) { return(NULL) }
      modules <- .res$res$ID
      mset    <- tmod_dbs[[.res$ds]][[.res$db]]$dbobj
      if(length(mset) < 2) {
        stop("Too few gene sets in the result list to show an upset plot")
      }
     #if(length(mset) > 50) {
     #  stop("Too many gene sets, use filter to make it smaller than 50")
     #}

      upset(modules, mset)
    }, width=800, height=800)

    output$tmodResTab <- renderDataTable({
      if(!isTruthy(res())) { return(NULL) }
      datatable(res()$res, escape=FALSE, selection='none', 
               options=list(pageLength=5, 
                            dom="Bfrtip", 
                            scrollX=TRUE)
                ) %>%
        formatSignif(columns=intersect(colnames(res()$res), 
                                       c("AUC", "cerno", "P.Value", "adj.P.Val")), digits=2)
    })

    output$table_sel_db <- renderUI({
      dbs <- names(tmod_res[[dataset()]][[1]])
      selectInput(NS(id, "db"), label="Database", choices=dbs, width="100%")
    })

    output$table_sel_sort <- renderUI({
      sorting <- names(tmod_res[[dataset()]][[1]][[1]])
      selectInput(NS(id, "sort"), label="Sorting", choices=sorting, width="100%")
    })

    observeEvent(input$select_button, {
      if(!is.null(gs_id)) {
        tmp <- unlist(strsplit(gsub("^go_", "", input$select_button), "-!-"))
        gs_id$ds <- tmp[1]
        gs_id$id <- tmp[2]
        gs_id$cntr <- tmp[3]
        gs_id$db <- tmp[4]
        gs_id$sort <- tmp[5]
      }
    })
  })
}

#' Launch a browser of tmod gene set enrichment analysis results
#'
#' Launch a shiny-based browser of tmod gene set enrichment analysis results
#'
#' To speed up launching the browser, you can load the tmod_dbs and
#' tmod_res objects first (using the functions get_tmod_dbs and
#' get_tmod_res).
#' @param pip pipeline object returned by `load_de_pipeline`
#' @param tmod_dbs tmod db object returned by get_tmod_dbs
#' @param tmod_res tmod results object returned by get_tmod_res
#' @param annot annotation data frame as returned by `get_annot`
#' @return does not return a value
#' @import dplyr
#' @importFrom shiny shinyApp renderText verbatimTextOutput textOutput renderUI uiOutput
#' @importFrom shiny tableOutput renderTable renderPlot plotOutput 
#' @importFrom shiny column fluidPage fluidRow mainPanel 
#' @importFrom shiny actionButton reactiveValues eventReactive
#' @importFrom shiny sidebarLayout sidebarPanel titlePanel tabPanel navbarPage updateNavbarPage tabsetPanel
#' @importFrom shiny selectInput numericInput sliderInput checkboxInput
#' @importFrom shiny downloadButton downloadHandler observeEvent reactiveVal isolate
#' @importFrom shiny showNotification removeNotification req numericInput
#' @importFrom shiny NS reactive is.reactive tagList moduleServer HTML h1 h2 h3 h4 br strong p
#' @importFrom shiny nearPoints hoverOpts brushedPoints
#' @importFrom shinyjs disable enable useShinyjs 
#' @importFrom grDevices dev.off pdf
#' @importFrom DT datatable formatSignif renderDataTable dataTableOutput
#' @importFrom colorDF summary_colorDF
#' @importFrom thematic thematic_shiny
#' @examples
#' \dontrun{
#' pip <- load_de_pipeline(config_file="DE_config.yaml")
#' tmod_dbs <- get_tmod_dbs(pip)
#' tmod_browser(pip, tmod_dbs)
#' }
#' @export
tmod_browser <- function(pip, tmod_dbs=NULL, tmod_res=NULL, annot=NULL) {

  message("preparing...")
  if(is.null(annot)) {
    message(" * Loading Annotation (consider using the annot option to speed this up)")
    annot  <- get_annot(pip)
  }

  if(is.null(tmod_res)) {
    message(" * Loading tmod results (consider using the tmod_res option to speed this up)")
    tmod_res  <- get_tmod_res(pip)
  }

  tmod_map <- get_tmod_mapping(pip)

  config <- get_config(pip)
  cntr_titles <- map_chr(config$contrasts$contrast_list, `[[`, "ID")
  names(cntr_titles) <- map_chr(config$contrasts$contrast_list, `[[`, "title")

  cntr   <- get_contrasts(pip) 

  if(is.null(tmod_dbs)) {
    message(" * Reading tmod_dbs (consider using it as an argument)...")
    tmod_dbs <- get_tmod_dbs(pip)
  }

  dbs <- names(tmod_dbs)
  sorting <- config$tmod$sort_by

  thematic_shiny(font="auto")

  ## prepare the tmod browser UI
  ui <- fluidPage(
    useShinyjs(),
    theme = bs_theme(bootswatch = "sandstone"),
    fluidRow(titlePanel(h1("tmod browser")), class="bg-primary"),
    fluidRow(HTML("<hr>")),
    tmodBrowserTableUI("tmod", cntr_titles),
    tmodBrowserPlotUI("tmodPlot"))
    
  server <- function(input, output, session) {
    gs_id <- reactiveValues()
    tmodBrowserTableServer("tmod", tmod_res, gs_id)
    tmodBrowserPlotServer("tmodPlot", gs_id, pip, tmod_dbs, tmod_map, cntr)
  }

  shinyApp(ui, server)
}


