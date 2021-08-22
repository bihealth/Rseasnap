## Construct the results table to display. Specifically, add action button
## for launching the plot.
.tmod_browser_prepare_res <- function(pip, but, tmod_res=NULL) {
  if(is.null(tmod_res)) {
    message(" * Reading tmod results (consider using it as an argument)...")
    tmod_res <- get_tmod_res(pip)
  }

  # prepare the tmod res
  tmod_res <- tmod_res %>% imap(~ {
    .cntr <- .y
    imap(.x, ~ {
           .dbname <- .y
           imap(.x, ~ {
                  .sort <- .y
                  .x %>% mutate(">"=sprintf(but, ID, .cntr, .dbname, .sort)) %>%
                    select(-cES, -cerno) %>%
                    arrange(P.Value) %>%
                    relocate(all_of(">"), .before=1)
            })
         })
  })

  return(tmod_res)
}


## create a datatable with the genes from a gene set
.tmod_browser_gene_table <- function(pip, but, id, db_name, cntr_name, sort_name, tmod_dbs, cntr, tmod_map) {

  message("Generating gene table")

  db <- tmod_dbs[[db_name]][["dbobj"]]
  genes <- db[["MODULES2GENES"]][[id]]
  genes_pid <- tmod_rev_db_map_ids(pip, ids=genes, dbname=db_name, tmod_dbs_mapping_obj=tmod_map)

  ret <- cntr[[cntr_name]] %>% filter(PrimaryID %in% genes_pid)
  ret <- ret %>% mutate('>' = sprintf(but, PrimaryID)) %>% relocate(all_of(">"), .before=1)
  message(sprintf("Gene table with %d genes", nrow(ret)))

  datatable(ret, escape=FALSE, selection='none',
                options=list(pageLength=5, dom="Bfrtip", scrollX=TRUE, buttons=c("copy", "csv", "excel"))) %>%
        formatSignif(columns=intersect(colnames(ret),
                                       c("baseMean", "log2FoldChange", "pvalue", "padj")), digits=2)
  
}


## given a module, contrast, sorting prepare a module info tab contents,
## including which genes are significant in the given contrast
.tmod_browser_mod_info <- function(pip, id, db_name, cntr_name, sort_name, tmod_dbs, cntr, tmod_map) {
  db <- tmod_dbs[[db_name]][["dbobj"]]
  ret <- sprintf("Module ID: %s\nDescription: %s\nContrast: %s\nDatabase: %s\nSorting: %s",
          id,
          db[["MODULES"]][id, ][["Title"]],
          cntr_name,
          db_name, 
          sort_name)

  genes <- db[["MODULES2GENES"]][[id]]
  genes_pid <- tmod_rev_db_map_ids(pip, ids=genes, dbname=db_name, tmod_dbs_mapping_obj=tmod_map)
  genes_sign <- genes_pid %in% (cntr[[cntr_name]] %>% filter(.data$padj < 0.05) %>% pull(.data$PrimaryID))
  ret <- paste0(ret, "\nSignificant genes:\n",     paste(genes[genes_sign],  collapse=","), "\n")
  ret <- paste0(ret, "\nNon-significant genes:\n", paste(genes[!genes_sign], collapse=","), "\n")
  return(ret)
}


#' @rdname tmodBrowserPlotServer
#' @export
tmodBrowserPlotUI <- function(id) {
    sidebarLayout(
      sidebarPanel(
        fluidRow(downloadButton(NS(id, "save"), "Save plot", class="bg-success")),
        fluidRow(verbatimTextOutput(NS(id, "modinfo"))),
        width=5
      ),
      mainPanel(
        fluidRow(verbatimTextOutput(NS(id, "cmdline"))),
        fluidRow(
                 tabsetPanel(
                             tabPanel("Plot", plotOutput(NS(id, "evidencePlot"))),
                             tabPanel("Genes", dataTableOutput(NS(id, "moduleGenes")))
                 )),
        width=7
      )
    )
}



## server module for viewing evidence plots

#' Shiny Module – tmod browser evidence plots
#'
#' Shiny Module – gene browser evidence plots
#' @param selmod a reactive
#' value (e.g. returned by tmodBrowserTableServer) and a list containing
#' the module id, tmod dataset id, contrast id and sorting type
#' @param pip pipeline object returned by `load_de_pipeline`
#' @param tmod_dbs tmod gene set databases returned by `get_tmod_dbs()`
#' @param tmod_map tmod gene set ID mapping returned by `get_tmod_mapping()`
#' @param id identifier (same as the one passed to geneBrowserTableUI)
#' @param cntr list of contrast results returned by `get_contrasts()`
#' @param annot data frame containing gene annotation 
#' @return returns a reactive value with a selected gene identifier
#' @export
tmodBrowserPlotServer <- function(id, selmod, pip, tmod_dbs, tmod_map, cntr, annot=NULL) {
  moduleServer(id, function(input, output, session) {
    message("Launching tmod plot server")

    ## gene_id is necessary to call the gene plots in another tab
    gene_id <- reactiveVal("")

    gene.but <- actionButton("go_%s", label=" \U25B6 ", 
                        onclick=sprintf('Shiny.onInputChange(\"%s-gene_select_button\",  this.id)', id),  
                        class = "btn-primary btn-sm")

    observeEvent(input$gene_select_button, {
      gene_id(gsub("^go_", "", input$gene_select_button))
    })

    disable("save")

    ## create the evidence plot and display the command line to replicate it
    output$evidencePlot <- renderPlot({
      mod <- req(selmod())
      enable("save")
      output$cmdline <- renderText({
        sprintf('plot_evidence(pip, id="%s", dbname="%s", contrast="%s")', 
                mod$id, mod$db, mod$cntr)
      })
      plot_evidence(pip, id=mod$id, dbname=mod$db, contrast=mod$cntr, tmod_dbs_obj=tmod_dbs)
    })


    output$modinfo <- renderText({
      mod <- req(selmod())
      ret <- .tmod_browser_mod_info(pip, 
                                    mod$id, mod$db, mod$cntr, mod$cntr, 
                                    tmod_dbs, cntr, tmod_map)
      return(ret)
    })

    output$moduleGenes <- renderDataTable({
      mod <- req(selmod())
      .tmod_browser_gene_table(pip, as.character(gene.but),
                                    mod$id, mod$db, mod$cntr, mod$cntr, 
                                    tmod_dbs, cntr, tmod_map)
    })
    


    ## save the plot as PDF
    output$save <- downloadHandler(
      filename = function() {
        mod <- req(selmod())
        ret <- sprintf("evidence_plot_%s_%s_%s.pdf", mod$db, mod$cntr, mod$id)
        return(ret)
      },
      content = function(file) {
        mod <- req(selmod())
        pdf(file=file, width=8, height=5)
        title <- sprintf("%s / %s\nContrast: %s / %s", 
                         mod$id, mod$db, mod$cntr, mod$sort)
        plot_evidence(pip, id=mod$id, dbname=mod$db, contrast=mod$cntr,
                      tmod_dbs_obj=tmod_dbs, tmod_dbs_mapping_obj=tmod_map, main=title)
        dev.off()
      }
    )
    gene_id
  })
}


#' @rdname tmodBrowserTableServer
#' @export
tmodBrowserTableUI <- function(id, cntr_titles, dbs, sorting) {

  but <- actionButton("uselessID", label=" \U25B6 ", class = "btn-primary btn-sm")

  ui <- sidebarLayout(
          sidebarPanel(
           fluidRow(selectInput(NS(id, "contrast"), label="Contrast", choices=cntr_titles, width="100%")),
           fluidRow(selectInput(NS(id, "db"),       label="Database", choices=dbs,         width="100%")),
           fluidRow(selectInput(NS(id, "sort"),     label="Sorting",  choices=sorting,     width="100%")),
           fluidRow(
             numericInput(NS(id, "f_auc"),  label="Filter by AUC", 
                          min=.5, max=1.0, step=0.1, value=0.65, width="50%"),
             numericInput(NS(id, "f_pval"), label="Filter by FDR", 
                          min=0, max=1.0, step=0.1, value=0.05, width="50%")
                    ),
           HTML(paste("Click on the", as.character(but), "buttons to view an evidence plot")),
           width=3
          ),
          mainPanel(
            dataTableOutput(NS(id, "tmodResTab")),
            width=9
          )
        )

  return(ui)
}

#' Shiny Module – tmod results browser table selection
#'
#' Shiny Module – tmod results browser table selection
#' @param tmod_res results of tmod analysis, returned by `get_tmod_res`
#' @param cntr_titles possibly named character vector with contrast names
#' @param dbs character vector with database names
#' @param sorting character vector with sorting options available
#' @param id identifier for the namespace of the module
#' @param pip pipeline object returned by `load_de_pipeline`
#' @return reactive value producing a list containing the module id, contrast id, db name and sort type.
#' @export
tmodBrowserTableServer <- function(id, pip, tmod_res) {


  moduleServer(id, function(input, output, session) {
    message("Launching tmod browser server")

    but <- actionButton("go_%s-!-%s-!-%s-!-%s", label=" \U25B6 ", 
                        onclick=sprintf('Shiny.onInputChange(\"%s-select_button\",  this.id)', id),  
                        class = "btn-primary btn-sm")

    tmod_res <- .tmod_browser_prepare_res(pip, as.character(but), tmod_res)

    output$tmodResTab <- renderDataTable({
      res <- tmod_res[[input$contrast]][[input$db]][[input$sort]] %>%
        filter(.data[["AUC"]] > input$f_auc & .data[["adj.P.Val"]] < input$f_pval)
      datatable(res, escape=FALSE, selection='none', options=list(pageLength=5)) %>%
        formatSignif(columns=intersect(colnames(res), 
                                       c("AUC", "cerno", "P.Value", "adj.P.Val")), digits=2)
    })


    ## this is the reactive expression returned by this module, known in
    ## tmodBrowserPlotServer as "selmod()"
    reactive({
      req(input$select_button)
      tmp <- unlist(strsplit(gsub("^go_", "", input$select_button), "-!-"))
      list(id=tmp[1], cntr=tmp[2], db=tmp[3], sort=tmp[4])
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
    annot  <- get_tmod_res(pip)
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
    tmodBrowserTableUI("tmod", cntr_titles, dbs, sorting),
    tmodBrowserPlotUI("tmodPlot"))
    
  server <- function(input, output, session) {
    selmod <- tmodBrowserTableServer("tmod", pip, tmod_res)
    tmodBrowserPlotServer("tmodPlot", selmod, pip, tmod_dbs, tmod_map, cntr)
  }

  shinyApp(ui, server)
}


