options(spinner.color="#47336F")
options(spinner.type=6)


## UI for the information about the pipeline
infoUI <- function(pipelines) {

  selector <- selectInput("select_pipeline",
                          label="Select pipeline:",
                          choices=pipelines)

  column(width=12,
    box(title=sprintf("%s: overview", textOutput("project_title")), width=6,
        solidHeader=TRUE, status="primary", collapsible=TRUE,
        fluidRow(selector),
        fluidRow(
           tableOutput("project_overview")
    )),
    box(title="Contrasts", width=6,
        solidHeader=TRUE, status="primary", collapsible=TRUE,
        fluidRow(
           tableOutput("contrasts_overview")
    )),
    box(title="Covariates:", width=12,
        solidHeader=TRUE, status="primary", collapsible=TRUE,
        column(width=11,
        fluidRow(
             dataTableOutput("covariates")
    )))
  )

}


## UI for the help system
helpUI <- function() {
  help_dir <- system.file("seapiper_manual/", package="Rseasnap")

  help_files <- list.files(help_dir, full.names = TRUE)

  ret <- lapply(help_files,
                function(f) {
                  title <- gsub("\\.md$", "", gsub("_", " ", basename(f)))
                  box(width=12, collapsible = TRUE, title = title,
                      solidHeader=TRUE, includeMarkdown(f), status="primary")
                })

  do.call(tagList, ret)
}


.pipeline_dashboard_header <- function(title) {

  dashboardHeader(title=img(src="icons/piper_horiz.png", alt="[seaPiper]"),
    tags$li(class="dropdown", h4(title, style="font-size:22px;color:white;padding-right:20px;"))
  )

}


.pipeline_dashboard_sidebar <- function() {
  dashboardSidebar(
       tags$head(
                 tags$link(rel = "stylesheet", type = "text/css", href = "css/seapiper.css"),
                 tags$link(rel="icon", type="image/png", sizes="32x32", href="icons/favicon-32x32.png"),
                 tags$link(rel="icon", type="image/png", sizes="16x16", href="icons/favicon-16x16.png"),
                 tags$link(rel="shortcut icon", type="image/x-icon", sizes="16x16", href="icons/favicon.ico")),
       sidebarMenu(
         # Setting id makes input$tabs give the tabName of currently-selected tab
         id = "navid",
         menuItem("Gene browser",  tabName = "gene_browser", icon = icon("dna")),
         menuItem("Tmod browser",  tabName = "tmod_browser", icon = icon("project-diagram")),
         menuItem("Disco plots",   tabName = "disco", icon = icon("chart-line")),
         menuItem("Panel plots",   tabName = "panel_plot", icon = icon("grip-vertical")),
         menuItem("PCA",           tabName = "pca", icon = icon("cube")),
         menuItem("Pipeline Info", tabName = "pip_info", icon = icon("info-circle")),
         menuItem("Help",          tabName = "help", icon = icon("question-circle"))
       )
  )

}


## prepare the actual tabs UI
.pipeline_dashboard_body <- function(data, title) {

  npip <- length(data[[1]])

  pipelines <- names(data[[1]])

  cntr_titles <- data[["cntr_titles"]]
  covar       <- data[["covar"]]

      t1 <- tabItem("gene_browser",
         box(title="Gene table", width=12, status="primary", 
             collapsible=TRUE,
             solidHeader=TRUE, geneBrowserTableUI("geneT", cntr_titles)),
         box(title="Gene info",  width=12, status="primary", 
             collapsible=TRUE,
             solidHeader=TRUE, geneBrowserPlotUI("geneP", contrasts=TRUE)),
        useShinyjs()
      )
    t2 <- tabItem("tmod_browser",
       box(title="Gene set enrichment overview", width=12, status="primary",
           collapsible=TRUE,
           solidHeader=TRUE, tmodBrowserTableUI("tmodT", cntr_titles)),
       box(title="Evidence plot", width=12, status="primary",
           collapsible=TRUE,
           solidHeader=TRUE, tmodBrowserPlotUI("tmodP"))
    )
#    t3 <- tabItem("disco",
#       box(title="Discordance / concordance plots", width=12, status="primary",
#       height="800px", solidHeader=TRUE, discoUI("disco", cntr_titles)),
#       )
#    t4 <- tabItem("panel_plot",
#       box(title="Panel plot", width=12, status="primary",
#           solidHeader=TRUE, tmodPanelPlotUI("panelP", dbs, sorting)))
#    t5 <- tabItem("pca",
#       box(title="Principal Component Analysis", width=12, status="primary",
#       solidHeader=TRUE, pcaUI("pca", covar, pca_names)),
#       useShinyjs()
#       )
      t6 <- tabItem("pip_info", 
         infoUI(pipelines))
         
      t7 <- tabItem("help", helpUI())
  dashboardBody(
    tabItems(
             t1, t2,# t3, t4, t5, 
             t6, t7
    ),
    style="min-height:1500px;"
  )
}


## load missing objects, prepare the data etc.
.prepare_data <- function(pip, annot=NULL, cntr=NULL, tmod_res=NULL, tmod_dbs=NULL, primary_id) {

  message("preparing...")
  if(is.null(annot))    { annot <- list() }
  if(is.null(cntr))     { cntr <- list() }
  if(is.null(tmod_res)) { tmod_res <- list() }
  if(is.null(tmod_dbs)) { tmod_dbs <- list() }

  data <- imap(pip, ~ {
    .pip <- .x
    .id  <- .y
    .prepare_data_single_pipeline(.id, .pip, annot, cntr, tmod_res, tmod_dbs, primary_id)
  })

  return(transpose(data))
}

## prepares the data structure for a single pipeline
.prepare_data_single_pipeline <- function(.id, .pip, annot, cntr, tmod_res, tmod_dbs, primary_id) {
  ret <- list()

  if(is.null(annot[[.id]])) {
    message(sprintf(" * Loading annotation for %s (consider using the annot option to speed this up)", .id))
    ret[["annot"]] <- get_annot(.pip)
  } else {
    ret[["annot"]] <- annot[[.id]]
  }

  if(is.null(cntr[[.id]])) {
    message(sprintf(" * Loading contrasts for %s (consider using the cntr option to speed this up)", .id))
    ret[["cntr"]] <- get_contrasts(.pip)
  } else {
    ret[["cntr"]] <- cntr[[.id]]
  }

  ret[["cntr"]] <- map(ret[["cntr"]], ~ .x %>% rownames_to_column(primary_id))

  if(is.null(tmod_res[[.id]])) {
    message(sprintf(" * Loading tmod results for %s (consider using the tmod_res option to speed this up)", .id))
    ret[["tmod_res"]] <- get_tmod_res(.pip)
  } else {
    ret[["tmod_res"]] <- tmod_res[[.id]]
  }

  if(is.null(tmod_dbs[[.id]])) {
    message(sprintf(" * Loading tmod databases for %s (consider using the tmod_dbs option to speed this up)", .id))
    ret[["tmod_dbs"]] <- get_tmod_dbs(.pip)
  } else {
    ret[["tmod_dbs"]] <- tmod_dbs[[.id]]
  }


  ret[["tmod_map"]] <- get_tmod_mapping(.pip)
  ret[["tmod_gl"]]  <- get_object(.pip, step="tmod", extension="gl.rds", as_list=TRUE)

  ret[["config"]]   <- get_config(.pip)
  ret[["covar"]]    <- get_covariates(.pip)

  ret[["annot_linkout"]] <- .prep_annot_linkout(ret[["annot"]], ret[["config"]])

  ret[["dbs"]]     <- names(tmod_dbs)
  ret[["sorting"]] <- config$tmod$sort_by

  ret[["rld"]]     <- get_object(.pip, step="DESeq2", extension="rld.blind.rds")
  ret[["rld"]]     <- ret[["rld"]]@assays@data@listData[[1]]

  ret[["pca"]] <- prcomp(t(ret[["rld"]]), scale.=TRUE)
  
  ret[["cntr_titles"]]        <- map_chr(ret[["config"]]$contrasts$contrast_list, `[[`, "ID")
  names(ret[["cntr_titles"]]) <- map_chr(ret[["config"]]$contrasts$contrast_list, `[[`, "title")
  ret[["cntr_titles"]]        <- ret[["cntr_titles"]][ ret[["cntr_titles"]] %in% names(ret[["cntr"]]) ]

  ret
}


#' Sea-snap pipeline browser
#'
#' Sea-snap pipeline browser
#' 
#' Launch a shiny app with gene and tmod results browsers for a pipeline.
#' Launching is much faster if the large objects (contrasts, tmod
#' databases) do not need to be loaded each time the pipeline browser is
#' started.
#' @param pip pipeline returned by `load_de_pipeline`
#' @param annot annotation returned by `get_annot`
#' @param cntr contrasts list returned by `get_contrasts`
#' @param tmod_dbs tmod databases returned by `get_tmod_dbs`
#' @param tmod_res tmod results returned by `get_tmod_res`
#' @param primary_id name of the column in the annotion data frame
#'        which corresponds to the primary gene identifier (including the
#'        row names of the contrasts results in the cntr object)
#' @param title Name of the pipeline to display
#' @importFrom shiny renderImage tags img icon imageOutput includeMarkdown
#' @importFrom shiny addResourcePath
#' @importFrom shinydashboard dashboardPage dashboardBody dashboardSidebar dashboardHeader
#' @importFrom shinydashboard box sidebarMenu tabItem tabItems updateTabItems menuItem
#' @examples
#' if(interactive()) {
#'   example_dir <- system.file("extdata/example_pipeline", package="Rseasnap")
#'   conf_f      <- file.path(example_dir, "DE_config.yaml")
#'   pip         <- load_de_pipeline(config_file = conf_f)
#'   pipeline_browser(pip)
#' }
#' @export
pipeline_browser <- function(pip, title="Pipeline browser", annot=NULL, cntr=NULL, tmod_res=NULL, tmod_dbs=NULL,
                             primary_id="PrimaryID") {

  addResourcePath("icons", system.file("icons", package="Rseasnap"))
  addResourcePath("css",   system.file("css", package="Rseasnap"))

  theme_set(theme_bw())

  ## pip can be a pipeline or a list of pipelines. In this first case, we
  ## change everything into a list.
  if(is(pip, "seasnap_DE_pipeline")) {
    pip <- list(ID1=pip)

    if(!is.null(annot))    { annot    <- list(ID1=annot)    }
    if(!is.null(cntr))     { cntr     <- list(ID1=cntr)     }
    if(!is.null(tmod_res)) { tmod_res <- list(ID1=tmod_res) }
    if(!is.null(tmod_dbs)) { tmod_dbs <- list(ID1=tmod_dbs) }
  }

  data <- .prepare_data(pip, annot, cntr, tmod_res, tmod_dbs, primary_id)

  thematic_shiny(font="auto")

  ## Prepare the UI
  header  <- .pipeline_dashboard_header(title)     
  sidebar <- .pipeline_dashboard_sidebar()
  body    <- .pipeline_dashboard_body(data, title)
  ui <- dashboardPage(header, sidebar, body, skin="purple", title=title)

  #   theme = bs_theme(primary = "#47336F", secondary = "#C6B3EB", 
  #                           font_scale = NULL, 
  #                           `enable-shadows` = TRUE, 
  #                           bootswatch = "united"),


  server <- function(input, output, session) {

    ## pipeline browser specific functions

#   output$project_overview   <- renderTable({ project_overview_table(config, title) })
#   output$contrasts_overview <- renderTable({ contrasts_overview_table(config) })
#   output$covariates         <- renderDataTable({ covariate_table(covar) })

    ## this reactive value holds the id of the selected gene, however the
    ## selection has been done
    gene_id <- reactiveVal()
    mod_id  <- reactiveVal()
    mod_id(list())

    gene_id1 <- geneBrowserTableServer("geneT", data[["cntr"]], data[["annot"]], 
                                       annot_linkout=data[["annot_linkout"]])
    gene_id3 <- tmodBrowserPlotServer("tmodP", mod_id, 
                                      tmod_dbs=data[["tmod_dbs"]], 
                                      cntr    =data[["cntr"]], 
                                      tmod_map=data[["tmod_map"]], 
                                      tmod_gl =data[["tmod_gl"]], 
                                      annot   =data[["annot"]])
#    gene_id2 <- discoServer("disco", cntr, annot)
    mod_id1  <- tmodBrowserTableServer("tmodT", data[["tmod_res"]], multilevel=TRUE)
#   mod_id2  <- tmodPanelPlotServer("panelP", cntr    =data[["cntr"]], 
#                                             tmod_res=data[["tmod_res"]],
#                                             tmod_dbs=data[["tmod_dbs"]], 
#                                             tmod_map=data[["tmod_map"]], 
#                                             annot   =data[["annot"]])

#    pcaServer("pca", pca$x, covar)
    
    ## combine events selecting a gene set 
    observeEvent(mod_id1(),  { 
                   mod_id(mod_id1()) 
    })
#    observeEvent(mod_id2(), { 
#      updateTabItems(session, "navid", "tmod_browser")
#      mod_id(mod_id2()) 
#    })

    ## combine events selecting a gene from gene browser and from disco
#    observeEvent(gene_id1(), { gene_id(gene_id1()) })
#    observeEvent(gene_id2(), { 
#      updateTabItems(session, "navid", "gene_browser")
#      gene_id(gene_id2())
#    })
#    observeEvent(gene_id3(), { 
#      updateTabItems(session, "navid", "gene_browser")
#      gene_id(gene_id3())
#    })

    geneBrowserPlotServer("geneP", gene_id1, covar=data[["covar"]], 
                          exprs=data[["rld"]], annot=data[["annot"]], cntr=data[["cntr"]])
  }

  shinyApp(ui, server)
}


## create a table with the overview of the project
project_overview_table <- function(config, title) {
  
  tmp1 <- data.frame(
                 c("Organism", "Taxon ID", "Formula", "Contrasts", "Tmod dbs", 
                   "Low count filter (absolute)", "Low count filter (samples)"),
                 c(config$organism$name, 
                   config$organism$taxon,
                   config$experiment$design_formula,
                   length(config$contrasts$contrast_list),
                   paste(map_chr(config$tmod$databases, `[[`, "title"), collapse=", "),
                   sprintf("Removed genes with < %d total counts", config$filter$low_counts),
                   sprintf("Kept genes with at least %d counts in at least %d samples",
                           config$filter$min_counts, config$filter$min_count_n)
                   ))
  colnames(tmp1) <- c("Title:", title)
  tmp1
}

contrasts_overview_table <- function(config) {

  data.frame(ID=map_chr(config$contrasts$contrast_list, `[[`, "ID"),
             Title=map_chr(config$contrasts$contrast_list, `[[`, "title")
             )


}

## create a datatable with covariates
covariate_table <- function(covar) {

  tmp <- covar[ , !colnames(covar) %in% c("filename", "md5") ]

  datatable(tmp, extensions="Buttons",
              options=list(dom="Bfrtip", scrollX=TRUE, 
                           buttons=c("copy", "csv", "excel")))

}


.prep_annot_linkout <- function(annot, config) {

  ret <- list()
  cn <- colnames(annot)

  if("ENSEMBL" %in% cn) {
    ret$ENSEMBL <- "https://www.ensembl.org/id/%s/"
  }
  if("ENSEMBLID" %in% cn) {
    ret$ENSEMBLID <- "https://www.ensembl.org/id/%s/"
  }
              
  if(config$organism$name == "human" & "SYMBOL" %in% cn) {
    ret$SYMBOL <- "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s"
  }

  if("ENTREZ" %in% cn) {
    ret$ENTREZ <- "https://www.ncbi.nlm.nih.gov/gene/?term=%s"
  }

  if("ENTREZID" %in% cn) {
    ret$ENTREZID <- "https://www.ncbi.nlm.nih.gov/gene/?term=%s"
  }

  if("REFSEQID" %in% cn) {
    ret$REFSEQID <- "https://www.ncbi.nlm.nih.gov/gene/?term=%s"
  }

  if("REFSEQ" %in% cn) {
    ret$REFSEQ <- "https://www.ncbi.nlm.nih.gov/gene/?term=%s"
  }

  ret
}
