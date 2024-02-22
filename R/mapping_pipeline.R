#' @rdname load_mapping_pipeline
#' @export
print.seasnap_mapping_pipeline <- function(x, ...) {

  cat("Sea-snap mapping pipeline object\n")
  cat(sprintf("Created on %s\n", x$timestamp))
}



#' Load a sea-snap mapping pipeline from a directory
#'
#' Load a sea-snap mapping pipeline from a directory
#'
#' This is the entry point for using a sea-snap mapping pipeline in R. First
#' create an object using this function, next use all the goodies that come
#' with Rseasnap.
#' 
#' There is only one way of creating the object: 
#' point to a mapping pipeline yaml configuration file.
#' @param x an object of class seasnap_mapping_pipeline
#' @param config_file path to the mapping yaml configuration file
#' @param ... Further arguments passed to the `print()` function
#' @return seasnap DE pipeline object of class "seasnap_DE_pipeline"
#' @importFrom yaml read_yaml
#' @importFrom lubridate now
#' @import colorDF
#' @importFrom utils read.table
#' @export
load_mapping_pipeline <- function(config_file="mapping_config.yaml") {

  ret <- list()
  class(ret) <- c("seasnap_mapping_pipeline", class(ret))
  ret$timestamp <- now()

  pip     <- read_yaml(config_file)
  ret$dir <- dirname(config_file)
  ret$pip_config  <- pip
  ret$config_file <- basename(config_file)

  #ret$out_path <- file.path(ret$dir, pip$pipeline_param$out_path_pattern)

  #ret$file_tab_path <- get_pipeline_path(ret$out_path, step="report", extension="tsv")
  #ret$file_tab      <- .get_table(ret$file_tab_path)
  #ret$file_tab$filename <- file.path(ret$dir, ret$file_tab$filename)

  #ret$config_path <- get_pipeline_path(ret$out_path, step="pipeline_report", extension="yaml")
  #ret$config      <- .get_yaml(ret$config_path)

  return(ret)
}


#' @rdname load_de_pipeline
#' @param object an object of class seasnap_DE_pipeline
#' @export
summary.seasnap_mapping_pipeline <- function(object, ...) {

  cat("Sea-snap mapping pipeline object\n")
  cat(sprintf(" created on %s\n", object$timestamp))
  #cat(sprintf(" output: %s\n", object$out_path))
  #cat(sprintf(" config: %s\n", object$config_path))
  #cat(sprintf(" file tab: %s\n", object$file_tab_path))
  #cat(sprintf(" design: %s\n", object$config$experiment$design_formula))
  #cat(sprintf("Pipeline file tab:\n"))
  #print(object$file_tab)


  # for(c in names(cntr)) {
  #   nsign <- sum(cntr[[c]]$padj < 0.05 & abs(cntr[[c]]$log2FoldChange) > 0.5, na.rm=TRUE)
  #   cat(sprintf(" %s: %d significant genes\n", c, nsign))
  # }

}


