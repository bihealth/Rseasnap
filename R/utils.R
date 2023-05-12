



## checks the pipeline object for validity
.check_de_obj <- function(x) {
  if(!"seasnap_DE_pipeline" %in% class(x)) {
    stop("x must be an object of class seasnap_DE_pipeline\n(use load_de_pipeline to create it)")
  }

  if(is.null(x$file_tab_path)) {
    stop("the object is corrupt, please use load_de_pipeline to create a new one")
  }

  if(file.exists(x$file_tab_path) && file.info(x$file_tab_path)$mtime[1] > attr(x$file_tab, "timestamp")) {
    warning("A newer file_tab file has been created\nConsider reloading with load_de_pipeline")
  }

}


## read and timestamp an RDS file
seasnap_readRDS <- function(x, path) {
  ts <- attr(x$config, "timestamp")
  if(!file.exists(path)) {
    stop(glue("File {path} does not exist"))
  }
  fi <- file.info(path)$mtime[1]

# if(ts < fi) {
#   warning(sprintf("Warning: file %s is newer than the pipeline config.\nYou should probably reload the pipeline with load_de_pipeline.", path))
# }

  ret <- readRDS(path)
  attr(ret, "timestamp") <- fi
  return(ret)
}

#' Create a seapiper_ds object
#'
#' Create an object that can be used with seaPiper
#'
#' @param x a seasnap_DE_pipeline object
#' @return An object of class seapiper_ds
#' @export
pip_to_seapiper <- function(x) {
  .check_de_obj(x)

  ret          <- list()
  ret$cntr     <- get_contrasts(x)
  ret$annot    <- get_annot(x)
  ret$tmod_dbs <- map(get_tmod_dbs(x), ~ .x$dbobj)
  ret$tmod_res <- get_tmod_res(x)
  ret$tmod_map <- get_tmod_mapping(x)
  ret$tmod_gl  <- get_tmod_gene_lists(x)

  ret$tmod_gl <- map(ret$tmod_gl, ~  # one for each of contrast
                             map(.[[1]], ~  # one for each sorting type
                                 match(names(.), ret$annot$PrimaryID) ))

  ret$pca      <- get_pca(x)
  ret$exprs    <- get_exprs(x)
  ret$config   <- get_config(x)
  ret$covar    <- get_covariates(x)
  #ret$cntr_titles <- get_contrast_names(x)
  ret[["cntr_titles"]]        <- map_chr(ret[["config"]]$contrasts$contrast_list, `[[`, "ID")
  names(ret[["cntr_titles"]]) <- map_chr(ret[["config"]]$contrasts$contrast_list, `[[`, "title")
  ret[["cntr_titles"]]        <- ret[["cntr_titles"]][ ret[["cntr_titles"]] %in% names(ret[["cntr"]]) ]

  message("returning seaPiper object")
  class(ret) <- c("seapiper_ds", class(ret))
  attr(ret, "primary_id") <- "PrimaryID"

  return(ret)
}

