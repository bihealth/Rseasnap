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


.tmod_rev_db_map_ids <- function(ids, dbname, tmod_dbs_mapping_obj) {

  mp_id <- tmod_dbs_mapping_obj$dbs[dbname]
  mapping <- tmod_dbs_mapping_obj$maps[[mp_id]]

  ret <- names(mapping)[ match(ids, mapping) ]
  
  if(all(is.na(ret))) {
    warning(glue("No IDs were found in the mapping... are you sure you are using the {dbname} IDs?\nCheck the annotation data frame (see `get_annot()`)"))
  }

  names(ret) <- ids

  return(ret)
}


## create a datatable with the genes from a gene set
.tmod_browser_gene_table <- function(but, ds, id, db_name, cntr_name, 
                                     sort_name, tmod_dbs, cntr, tmod_map, primary_id) {

  db <- tmod_dbs[[db_name]][["dbobj"]]
  genes <- db[["MODULES2GENES"]][[id]]
  genes_pid <- .tmod_rev_db_map_ids(ids=genes, dbname=db_name, tmod_dbs_mapping_obj=tmod_map)

  ret <- cntr[[cntr_name]] %>% filter(.data[[primary_id]] %in% genes_pid)
  ret <- ret %>% mutate('>' = sprintf(but, ds, .data[[primary_id]])) %>% relocate(all_of(">"), .before=1) %>%
    arrange(.data[["pvalue"]])

  num_cols <- colnames(ret)[ map_lgl(ret, is.numeric) ]

  datatable(ret, escape=FALSE, selection='none',
                options=list(pageLength=5, dom="Bfrtip", scrollX=TRUE, buttons=c("copy", "csv", "excel"))) %>%
        formatSignif(columns=num_cols, digits=2)
  
}

.plot_evidence <- function(mod_id, cntr_id, db_id, sort_id, cntr, tmod_dbs,
                           tmod_gl=NULL, tmod_map=NULL, annot=NULL, primary_id) {

  mset <- tmod_dbs[[db_id]][["dbobj"]]
  cntr <- cntr[[ cntr_id ]]

  if(!all(c("pvalue", "log2FoldChange", "padj") %in% colnames(cntr))) {
    stop("Parameter `cntr` must have columns log2FoldChange, pvalue and padj")
  }

  if(!primary_id %in% colnames(cntr)) {
    cntr[[primary_id]] <- rownames(cntr)
  }

  if(is.null(tmod_gl)) {
    ## create an ad hoc list
    cntr <- cntr[ order(cntr$pvalue), ]
    gl <- tmod_map$maps[[ tmod_map$dbs[[ db_id ]] ]][ cntr[[primary_id]] ]
  } else {
    gl <- tmod_gl[[cntr_id]][[db_id]][[sort_id]]
  }

  symbols <- names(gl)

  if("SYMBOL" %in% colnames(annot)) {
    symbols <- annot[["SYMBOL"]][ match(symbols, annot[[primary_id]]) ]
  }

  names(symbols) <- names(gl)

  ## We need to figure out what genes are in the gene set to select them as
  ## labels

  genes <- getModuleMembers(mod_id, mset)[[mod_id]]
  if(is.null(genes)) {
    stop(sprintf("Gene set %s not found in db %s", mod_id, db_id))
  }

  sel                <- gl %in% genes
  gene.labels        <- symbols[sel]
  names(gene.labels) <- gl[sel]
  
  ## now for some color
  # make sure that cntr is in the same order as gl
  cntr <- cntr[ match(names(gl), cntr[[primary_id]]), ]

  p <- cntr$padj
  l <- cntr$log2FoldChange

  colors <- ifelse(l < 0, 
      ifelse(p < .05, 'blue', '#000066'),
      ifelse(p < .05, 'red', '#660000'))
  names(colors) <- gl

  evidencePlot(gl, m=mod_id, mset=mset, gene.colors=colors, gene.labels=gene.labels)
}

## given a module, contrast, sorting: prepare a module info tab contents,
## including which genes are significant in the given contrast
.tmod_browser_mod_info <- function(id, ds, db_name, cntr_name, sort_name, tmod_dbs, cntr, tmod_map) {
  db <- tmod_dbs[[db_name]][["dbobj"]]
  ret <- sprintf("Module ID: %s\nDescription: %s\nData set: %s\nContrast: %s\nDatabase: %s\nSorting: %s",
          id,
          db[["MODULES"]][id, ][["Title"]],
          ds,
          cntr_name,
          db_name, 
          sort_name)

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
                             tabPanel("Plot", withSpinner(plotOutput(NS(id, "evidencePlot"), height="100%"))),
                             tabPanel("Genes", withSpinner(dataTableOutput(NS(id, "moduleGenes"))))
                 )),
        width=7
      )
    )
}



## server module for viewing evidence plots

#' Shiny Module – tmod browser evidence plots
#'
#' Shiny Module – gene browser evidence plots
#'
#' This part is a bit complex, because a lot of different data go into
#' the evidence plots. To create a plot, following elements are necessary:
#'
#'  * ordered gene list: this is the same gene list that is used as an
#'    input to tmod
#'  * a list of tmod gene set databases for which the enrichments have been run
#'  * contrast data frame – for displaying gene names in color
#'  * if no gene list is given, then using the contrast data frame it is
#'    possible to create a list ordered by p-values. However, since the
#'    gene set database might use a different type of identifiers than the
#'    PrimaryID column of the contrasts data frame, it is necessary to
#'    provide a mapping between the PrimaryIDs and the database gene IDs as well.
#'
#' 
#'
#' **tmod gene lists**. To display an evidence plot, we need to have an
#' ordered list of genes. This has to be provided from outside, as many
#' different sortings are possible. The parameter tmod_gl is a multi level 
#' list:
#'  * First level: contrasts
#'  * Second level: gene set databases
#'  * Third level: sorting type
#' 
#' For example, then the `tmod_gl[["Contrast1"]][["tmod"]][["pval"]]` is a
#' vector of gene identifiers which come from the tmod database (in this
#' case, gene names). The vector is named and the names are primary IDs
#' matching the contrast data frames.
#'
#' If this argument is NULL, then the genes will be ordered by p-values
#' from the contrast object provided. However, in this case it is necessary
#' to provide a mapping between the PrimaryIDs of the contrasts and the
#' gene identifiers used by the gene set database.
#'
#' @param id ID of the shiny module
#' @param selmod a reactive
#' value (e.g. returned by tmodBrowserTableServer) and a list containing
#' the module id, tmod dataset id, contrast id and sorting type
#' @param tmod_dbs tmod gene set databases returned by `get_tmod_dbs()`
#' @param tmod_map tmod gene set ID mapping returned by `get_tmod_mapping()`
#' @param tmod_gl tmod gene lists. See details.
#' @param id identifier (same as the one passed to geneBrowserTableUI)
#' @param primary_id name of the column which holds the primary identifiers
#' @param cntr list of contrast results returned by `get_contrasts()`
#' @param annot data frame containing gene annotation 
#' @return returns a reactive value with a selected gene identifier
#' @importFrom shinyBS popify
#' @export
tmodBrowserPlotServer <- function(id, selmod, tmod_dbs, cntr, tmod_map=NULL, tmod_gl=NULL, annot=NULL, 
                                  primary_id="PrimaryID") {

  stopifnot(!is.null(tmod_gl) || !is.null(tmod_map))

  if(!is.data.frame(annot)) {
    message("tmodBrowserPlotServer: running in multilevel mode")
  } else {
    tmod_dbs <- list(default=tmod_dbs)
    cntr     <- list(default=cntr)
    tmod_map <- list(default=tmod_map)
    tmod_gl  <- list(default=tmod_gl)
    annot    <- list(default=annot)
  }
    
  moduleServer(id, function(input, output, session) {
    message("Launching tmod plot server")

    ## gene_id is necessary to call the gene plots in another tab
    gene_id <- reactiveVal("")

    gene.but <- actionButton("go~%s~%s", label=" \U25B6 ", 
                        onclick=sprintf('Shiny.onInputChange(\"%s-gene_select_button\",  this.id)', id),  
                        class = "btn-primary btn-sm")

    gene_id <- reactiveVal(list(ds=NA, id=NA))

    observeEvent(input$gene_select_button, {
      .ids <- strsplit(input$gene_select_button, '~')[[1]]
      gene_id(list(ds=.ids[2], id=.ids[3]))
    })

    disable("save")

    ## create the evidence plot and display the command line to replicate it
    output$evidencePlot <- renderPlot({
      message("Rendering plot")
      mod <- req(selmod())

      if(is.null(mod$id)) { return(NULL) }
      enable("save")

      ds <- mod$ds
      .plot_evidence(mod_id=mod$id, cntr_id=mod$cntr, db_id=mod$db, sort_id=mod$sort, 
                     cntr=cntr[[ds]], tmod_dbs=tmod_dbs[[ds]], tmod_gl=tmod_gl[[ds]], 
                     tmod_map=tmod_map[[ds]], annot=annot[[ds]], primary_id)
    }, width=800, height=600)


    output$modinfo <- renderText({
      mod <- req(selmod())
      if(is.null(mod$id)) { return(NULL) }
      ret <- .tmod_browser_mod_info(mod$id, mod$ds, mod$db, mod$cntr, mod$sort, 
                                    tmod_dbs[[mod$ds]], cntr[[mod$ds]], tmod_map[[mod$ds]])
      return(ret)
    })

    output$moduleGenes <- renderDataTable({
      mod <- req(selmod())
      ds <- mod$ds
      .tmod_browser_gene_table(as.character(gene.but), ds,
                                    mod$id, mod$db, mod$cntr, mod$cntr, 
                                    tmod_dbs[[ds]], cntr[[ds]], tmod_map[[ds]], primary_id)
    })
    


    ## save the plot as PDF
    output$save <- downloadHandler(
      filename = function() {
        mod <- req(selmod())
        if(is.null(mod$id)) { return(NULL) }
        ret <- sprintf("evidence_plot_%s_%s_%s_%s.pdf", mod$ds, mod$db, mod$cntr, mod$id)
        return(ret)
      },
      content = function(file) {
        mod <- req(selmod())
        if(is.null(mod$id)) { return(NULL) }
        pdf(file=file, width=8, height=5)
        title <- sprintf("%s / %s\nContrast: %s / %s", 
                         mod$id, mod$db, mod$cntr, mod$sort)

        ds <- mod$ds
        .plot_evidence(mod_id=mod$id, cntr_id=mod$cntr, db_id=mod$db, sort_id=mod$sort, 
                     cntr=cntr[[ds]], tmod_dbs=tmod_dbs[[ds]], tmod_gl=tmod_gl[[ds]], 
                     tmod_map=tmod_map[[ds]], annot=annot[[ds]], primary_id)
        dev.off()
      }
    )
    gene_id
  })
}


#' @rdname tmodBrowserTableServer
#' @export
tmodBrowserTableUI <- function(id, cntr_titles) {

  cntr_titles <- .prep_cntr_titles(cntr_titles)

  but <- actionButton("uselessID", label=" \U25B6 ", class = "btn-primary btn-sm")

  ui <- sidebarLayout(
          sidebarPanel(
           fluidRow(selectInput(NS(id, "contrast"), label="Contrast", choices=cntr_titles, width="100%")),
           fluidRow(uiOutput(NS(id, "table_sel_db"))),
           fluidRow(uiOutput(NS(id, "table_sel_sort"))),
#           fluidRow(selectInput(NS(id, "db"),       label="Database", choices=dbs,         width="100%")),
#           fluidRow(selectInput(NS(id, "sort"),     label="Sorting",  choices=sorting,     width="100%")),
           fluidRow(
             checkboxInput(NS(id, "filter"), label="Filter results", value=TRUE),
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
#' @param id identifier for the namespace of the module
#' @param multilevel if TRUE, the results are grouped in data sets
#' @return reactive value producing a list containing the module id, contrast id, db name and sort type.
#' @export
tmodBrowserTableServer <- function(id, tmod_res, multilevel=FALSE) {

  if(!multilevel) {
    tmod_res <- list(default=tmod_res)
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

    output$tmodResTab <- renderDataTable({
      if(!(
           isTruthy(dataset()) &&
           isTruthy(contrast()) &&
           isTruthy(input$db) &&
           isTruthy(input$sort)
         )) { return(NULL) }
         
      res <- tmod_res[[dataset()]][[contrast()]][[input$db]][[input$sort]] 
      
      if(input$filter) {
        res <- res %>%
        filter(.data[["AUC"]] > input$f_auc & .data[["adj.P.Val"]] < input$f_pval)
      }

      datatable(res, escape=FALSE, selection='none', options=list(pageLength=5)) %>%
        formatSignif(columns=intersect(colnames(res), 
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

    ## this is the reactive expression returned by this module, known in
    ## tmodBrowserPlotServer as "selmod()"
    selmod <- reactiveVal()

    observeEvent(input$select_button, {
      tmp <- unlist(strsplit(gsub("^go_", "", input$select_button), "-!-"))
      selmod(list(ds=tmp[1], id=tmp[2], cntr=tmp[3], db=tmp[4], sort=tmp[5]))
    })

    return(selmod)

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
    selmod <- tmodBrowserTableServer("tmod", tmod_res)
    tmodBrowserPlotServer("tmodPlot", selmod, pip, tmod_dbs, tmod_map, cntr)
  }

  shinyApp(ui, server)
}


