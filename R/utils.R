



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



