


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
#' @examples
#' if(interactive) {
#'   example_dir <- system.file("extdata/example_pipeline", package="Rseasnap")
#'   conf_f      <- file.path(example_dir, "DE_config.yaml")
#'   pip         <- load_de_pipeline(config_file = conf_f)
#'   pipeline_browser(pip)
#' }
#' @export
pipeline_browser <- function(pip, annot=NULL, cntr=NULL, tmod_res=NULL, tmod_dbs=NULL,
                             primary_id="PrimaryID") {
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

  ui <- navbarPage(
      "Sea-snap pipeline browser",
      id = "navid",
      theme = bs_theme(primary = "#47336F", secondary = "#C6B3EB", 
                              font_scale = NULL, 
                              `enable-shadows` = TRUE, 
                              bootswatch = "united"),
        tabPanel(h3("Gene browser"),
            fluidRow(verbatimTextOutput("msg")),
            geneBrowserTableUI("geneT", cntr_titles),
            geneBrowserPlotUI("geneP", covar, contrasts=TRUE),
            value="gbrowser"
          ),
        tabPanel(h3("tmod browser"),
            tmodBrowserTableUI("tmodT", cntr_titles, dbs, sorting),
            tmodBrowserPlotUI("tmodP"),
            value="tbrowser"
        ),
        tabPanel(h3("disco"),
            discoUI("disco", cntr_titles),
            value="disco"
        ),
        tabPanel(h3("PCA"),
            pcaUI("pca", covar, colnames(pca$x)),
            value="pca"
        ),
        useShinyjs()
    )

  server <- function(input, output, session) {
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
      updateNavbarPage(session, "navid", "gbrowser")
      gene_id(gene_id2())
    })

    geneBrowserPlotServer("geneP", gene_id, covar=covar, 
                          exprs=rld, annot=annot, cntr=cntr)
  }

  shinyApp(ui, server)
}


