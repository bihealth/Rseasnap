get_pipeline_path <- function(path_pattern, step, extension, contrast="all") {
  path <- gsub("\\{extension\\}", extension, gsub("\\{step\\}", step, gsub("\\{contrast\\}", contrast, path_pattern)))
  return(path)
}



#' @rdname load_de_pipeline
#' @export
print.seasnap_DE_pipeline <- function(x, ...) {

  cat("Sea-snap DE pipeline object\n")
  cat(sprintf("Created on %s\n", x$timestamp))
}



#' Load a sea-snap DE pipeline from a directory
#'
#' Load a sea-snap DE pipeline from a directory
#'
#' This is the entry point for using a sea-snap DE pipeline in R. First
#' create an object using this function, next use all the goodies that come
#' with Rseasnap.
#' 
#' There are two ways of creating the object: either use an existing
#' configuration object
#' (config object in R snippets) or point to a DE pipeline yaml configuration file.
#' In any case, the config is only consulted to get the output file path of
#' the pipeline; then, the files file_tab and config are read to ensure that
#' the object is up to date.
#' @param path path to the pipeline directory. This is overridden if config
#'        and file_tab are defined.
#' @param config object accessible in R snippets of the DE
#'         pipeline containing the pipeline configuration and produced files. If tis
#'         parameter is defined, it overrides the `path` parameter.
#' @return seasnap DE pipeline object of class "seasnap_DE_pipeline"
#' @importFrom yaml read_yaml
#' @import lubridate
#' @import colorDF
#' @export
load_de_pipeline <- function(config_file="DE_config.yaml", config=NULL) {

  ret <- list()
  class(ret) <- c(class(ret), "seasnap_DE_pipeline")
  ret$timestamp <- now()

  if(!is.null(config)) {
    if(is.null(config$pipeline_param) || is.null(config$pipeline_param$out_path_pattern)) {
      stop("config does not contain config$pipeline_param$out_path_pattern")
    }
    pip <- config
  } else {
    pip <- read_yaml(config_file)
  }

  ret$out_path <- pip$pipeline_param$out_path_pattern

  ret$file_tab_path <- get_pipeline_path(ret$out_path, step="report", extension="tsv")
  ret$file_tab      <- .get_table(ret$file_tab_path)

  ret$config_path <- get_pipeline_path(ret$out_path, step="pipeline_report", extension="yaml")
  ret$config      <- .get_yaml(ret$config_path)

  return(ret)
}


## load a file and make a timestamp
.get_file_timestamped <- function(path, func, oldobj=NULL, check_for_newer=TRUE, ...) {

  if(!file.exists(path)) {
    stop(glue("File {path} does not exist."))
  }

  timestamp <- file.info(path)$mtime[1]

  if(!is.null(oldobj) && !is.null(attr(oldobj, "timestamp")) && check_for_newer) {
    if(!attr(oldobj, "timestamp") < timestamp) {
      warning("re-using old object")
      return(oldobj)
    }
  }

  ret <- func(path, ...)
  attr(ret, "timestamp") <- timestamp
  return(ret)
}

## load a data frame
.get_table_func <- function(path) {
  ret <- as.colorDF(read.table(path, header=TRUE, stringsAsFactors=FALSE))
  return(ret)
}

## load a data frame, timestamped
.get_table <- function(path, oldobj=NULL, check_for_newer=TRUE) {
  .get_file_timestamped(path, .get_table_func, oldobj, check_for_newer)
}

## load a yaml file, timestamped
.get_yaml <- function(path, oldobj=NULL, check_for_newer=TRUE) {
  .get_file_timestamped(path, yaml.load_file, oldobj, check_for_newer)
}




