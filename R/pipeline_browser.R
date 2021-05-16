


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
#' @export
pipeline_browser <- function(pip, annot=NULL, cntr=NULL, tmod_res=NULL, tmod_dbs=NULL) {

  message("preparing...")
  if(is.null(annot)) {
    message(" * Loading Annotation (consider using the annot option to speed this up)")
    annot  <- get_annot(pip)
  }

  # prepare the contrast tables
  if(is.null(cntr)) {
    message(" * Loading contrasts... (consider using the cntr option to speed this up)")
    cntr <- get_contrasts(pip)
  }

  if(is.null(tmod_res)) {
    message(" * Loading tmod results (consider using the tmod_res option to speed this up)")
    tmod_res  <- get_tmod_res(pip)
  }

  tmod_map <- get_tmod_mapping(pip)
  config   <- get_config(pip)
  covar    <- get_covariates(pip)

  if(is.null(tmod_dbs)) {
    message(" * Reading tmod_dbs (consider using it as an argument)...")
    tmod_dbs <- get_tmod_dbs(pip)
  }

  dbs <- names(tmod_dbs)
  sorting <- config$tmod$sort_by

  all_covars    <- covar %>% summary_colorDF() %>% filter(unique > 1) %>% pull(.data$Col)
  default_covar <- .default_covar(covar, all_covars, default="group")

  rld    <- get_object(pip, step="DESeq2", extension="rld.blind.rds")
  rld    <- rld@assays@data@listData[[1]]
  
  cntr_titles <- map_chr(config$contrasts$contrast_list, `[[`, "ID")
  names(cntr_titles) <- map_chr(config$contrasts$contrast_list, `[[`, "title")

  thematic_shiny(font="auto")

  ui <- navbarPage(
      "Sea-snap pipeline browser",
      theme = bs_theme(bootswatch = "united"),
        tabPanel(h3("Gene browser"),
            geneBrowserTableUI("geneT", cntr_titles),
            geneBrowserPlotUI("geneP", covar)
          ),
        tabPanel(h3("tmod browser"),
            tmodBrowserTableUI("tmodT", cntr_titles, dbs, sorting),
            tmodBrowserPlotUI("tmodP")
        ),
        useShinyjs()
    )
  server <- function(input, output, session) {
    gene_id <- geneBrowserTableServer("geneT", cntr, annot)
    geneBrowserPlotServer("geneP", gene_id, covar, rld, annot)
    mod_id <- tmodBrowserTableServer("tmodT", pip, tmod_res)
    tmodBrowserPlotServer("tmodP", mod_id, pip, tmod_dbs, tmod_map, cntr)
  }

  shinyApp(ui, server)
}
