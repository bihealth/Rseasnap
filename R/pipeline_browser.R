infoUI <- function(title) {

         #box(title="Principal Component Analysis", width=12, status="primary",
         #solidHeader=TRUE, pcaUI("pca", covar, colnames(pca$x))),
         #useShinyjs()
  

  column(width=12,
    box(title=sprintf("%s: overview", title), width=6,
        solidHeader=TRUE, status="primary", collapsible=TRUE,
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

  message("preparing...")
  if(is.null(annot)) {
    message(" * Loading Annotation (consider using the annot option to speed this up)")
    annot  <- get_annot(pip)
  }

  # prepare the contrast tables
  if(is.null(cntr)) {
    message(" * Loading contrasts (consider using the cntr option to speed this up)")
    cntr <- get_contrasts(pip)
  }

  cntr <- map(cntr, ~ .x %>% rownames_to_column(primary_id))


  if(is.null(tmod_res)) {
    message(" * Loading tmod results (consider using the tmod_res option to speed this up)")
    tmod_res  <- get_tmod_res(pip)
  }

  tmod_map <- get_tmod_mapping(pip)
  config   <- get_config(pip)
  covar    <- get_covariates(pip)

  if(is.null(tmod_dbs)) {
    message(" * Loading tmod_dbs (consider using the tmod_rdbses option to speed this up)")
    tmod_dbs <- get_tmod_dbs(pip)
  }

  dbs <- names(tmod_dbs)
  sorting <- config$tmod$sort_by

  rld    <- get_object(pip, step="DESeq2", extension="rld.blind.rds")
  rld    <- rld@assays@data@listData[[1]]

  pca <- prcomp(t(rld), scale.=TRUE)
  
  cntr_titles <- map_chr(config$contrasts$contrast_list, `[[`, "ID")
  names(cntr_titles) <- map_chr(config$contrasts$contrast_list, `[[`, "title")
  cntr_titles <- cntr_titles[ cntr_titles %in% names(cntr) ]

  thematic_shiny(font="auto")

  header <- dashboardHeader(title=img(src="icons/piper_horiz.png", alt="[seaPiper]"),
    tags$li(class="dropdown", h4(title, style="font-size:22px;color:white;padding-right:20px;"))
  )
     
  sidebar <- dashboardSidebar(
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
         menuItem("PCA",           tabName = "pca", icon = icon("cube")),
         menuItem("Pipeline Info", tabName = "pip_info", icon = icon("info-circle")),
         menuItem("Help",          tabName = "help", icon = icon("question-circle"))
       )
     )
     
  #   theme = bs_theme(primary = "#47336F", secondary = "#C6B3EB", 
  #                           font_scale = NULL, 
  #                           `enable-shadows` = TRUE, 
  #                           bootswatch = "united"),

  body <- dashboardBody(
    tabItems(
      tabItem("gene_browser",
         box(title="Gene table", width=12, status="primary", 
             collapsible=TRUE,
             solidHeader=TRUE, geneBrowserTableUI("geneT", cntr_titles)),
         box(title="Gene info",  width=12, status="primary", 
             collapsible=TRUE,
             solidHeader=TRUE, geneBrowserPlotUI("geneP", covar, contrasts=TRUE))
      ),
      tabItem("tmod_browser",
         box(title="Gene set enrichment overview", width=12, status="primary",
             collapsible=TRUE,
             solidHeader=TRUE, tmodBrowserTableUI("tmodT", cntr_titles, dbs, sorting)),
         box(title="Evidence plot", width=12, status="primary",
             collapsible=TRUE,
             solidHeader=TRUE, tmodBrowserPlotUI("tmodP"))
      ),
      tabItem("disco",
         box(title="Discordance / concordance plots", width=12, status="primary",
         height="800px", solidHeader=TRUE, discoUI("disco", cntr_titles)),
         ),
      tabItem("pca",
         box(title="Principal Component Analysis", width=12, status="primary",
         solidHeader=TRUE, pcaUI("pca", covar, colnames(pca$x))),
         useShinyjs()
         ),
      tabItem("pip_info", 
         infoUI(title)),
         
      tabItem("help",
         helpUI())
    )
  )

  ui <- dashboardPage(header, sidebar, body, skin="purple")

  server <- function(input, output, session) {

    ## pipeline browser specific functions

    output$project_overview <- renderTable({ project_overview_table(config, title) })
    output$contrasts_overview <- renderTable({ contrasts_overview_table(config) })
    output$covariates <- renderDataTable({ covariate_table(covar) })

    ## this reactive value holds the id of the selected gene, however the
    ## selection has been done
    gene_id <- reactiveVal()

    gene_id1 <- geneBrowserTableServer("geneT", cntr, annot)
    mod_id   <- tmodBrowserTableServer("tmodT", pip, tmod_res)
    tmodBrowserPlotServer("tmodP", mod_id, pip, tmod_dbs, tmod_map, cntr)
    gene_id2 <- discoServer("disco", cntr, annot)

    pcaServer("pca", pca$x, covar)

    ## combine events selecting a gene from gene browser and from disco
    observeEvent(gene_id1(), { gene_id(gene_id1()) })
    observeEvent(gene_id2(), { 
      updateTabItems(session, "navid", "gene_browser")
      gene_id(gene_id2())
    })

    geneBrowserPlotServer("geneP", gene_id, covar=covar, 
                          exprs=rld, annot=annot, cntr=cntr)
  }

  shinyApp(ui, server)
}


## create a table with the overview of the project
project_overview_table <- function(config, title) {
  
  tmp1 <- data.frame(
                 c("Organism", "Taxon ID", "Formula", "Contrasts", "Tmod dbs"),
                 c(config$organism$name, 
                   config$organism$taxon,
                   config$experiment$design_formula,
                   length(config$contrasts$contrast_list),
                   paste(map_chr(config$tmod$databases, `[[`, "title"), collapse=", ")
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

