get_pipeline_path <- function(path_pattern, step, extension, contrast="all") {
  path <- gsub("\\{extension\\}", extension, gsub("\\{step\\}", step, gsub("\\{contrast\\}", contrast, path_pattern)))
  return(path)
}




print.seasnap_DE_pipeline <- function(x, ...) {

  cat("Sea-snap DE pipeline object\n")
  cat(sprintf("Created on %s\n", x$timestamp))
}



#' Load a sea-snap DE pipeline from a directory
#'
#' Load a sea-snap DE pipeline from a directory
#'
#' @param path path to the pipeline directory
#' @return seasnap DE pipeline object of class "seasnap_DE_pipeline"
#' @importFrom yaml read_yaml
#' @import lubridate
#' @import colorDF
#' @export
load_de_pipeline <- function(config_file="DE_config.yaml") {

  ret <- list()
  class(ret) <- c(class(ret), "seasnap_DE_pipeline")

  ret$timestamp <- now()

  pip <- read_yaml(config_file)

  ret$out_path <- pip$pipeline_param$out_path_pattern
  ret$log_path <- pip$pipeline_param$log_path_pattern

  ret$file_tab_path <- get_pipeline_path(ret$out_path, step="report", extension="tsv")
  ret$file_tab      <- .get_table(ret$file_tab_path)

  ret$config_path <- get_pipeline_path(ret$out_path, step="pipeline_report", extension="yaml")
  ret$config      <- .get_yaml(ret$config_path)

  return(ret)
}


.get_file_timestamped <- function(path, func, oldobj=NULL, check_for_newer=TRUE, ...) {

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


.get_table_func <- function(path) {
  ret <- as.colorDF(read.table(path, header=TRUE, stringsAsFactors=FALSE))
  return(ret)
}

.get_table <- function(path, oldobj=NULL, check_for_newer=TRUE) {
  .get_file_timestamped(path, .get_table_func, oldobj, check_for_newer)
}

.get_yaml <- function(path, oldobj=NULL, check_for_newer=TRUE) {
  .get_file_timestamped(path, yaml.load_file, oldobj, check_for_newer)
}




