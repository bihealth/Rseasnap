#' Get the file tab
#'
#' Get a tabular representation of all output files created by the pipeline
#' @seealso load_de_pipeline
#' @param x an object of class seasnap_DE_pipeline
#' @return A color data frame with files
#' @export
get_file_tab <- function(x) {
  .check_de_obj(x)

  return(x$file_tab)
}

#' Get all contrast names
#'
#' Get the names of all defined contrasts from a pipeline
#' @param x an object of class seasnap_DE_pipeline
#' @return A character vector with contrast names
#' @export
get_contrast_names <- function(x) {
  .check_de_obj(x)

  ret <- setdiff(unique(x$file_tab$contrast), "all")
  ord <- order(gsub(".*_(ID[0-9]*)$", "\\1", ret))
  return(ret[ord])
}


#' Get all results of DE analysis
#'
#' Get results of the DE analysis for all contrasts (default) or a
#' selection of the contrasts.
#' @param x an object of class seasnap_DE_pipeline
#' @param contrasts an optional character vector of the contrasts to return
#' @return a named list; each object of the list is a data frame with the results of DE analysis of a given contrast.
#' @export
get_contrasts <- function(x, contrasts=NULL) {
  .check_de_obj(x)

  if(is.null(contrasts)) {
    contrasts <- get_contrast_names(x)
  }

  ret <- lapply(contrasts, function(c) {
    sel <- with(x$file_tab, step == "contrast" & contrast == c & extension == "rds")
    path <- x$file_tab$filename[sel]
    as.data.frame(seasnap_readRDS(x, path))
  })

  names(ret) <- contrasts
  return(ret)
}

#' Get all counts from the pipeline
#'
#' Get counts from the import_gene_counts step
#' @param x an object of class seasnap_DE_pipeline
#' @param as_matrix if FALSE, a data frame will be returned
#' @return either a matrix with raw counts, rownames corresponding to PrimaryIDs
#'         and column names corresponding to sample names, or data frame
#'         with first column being the PrimaryIDs
#' @importFrom utils read.csv
#' @export
get_counts <- function(x, as_matrix=TRUE) {
  cfile <- get_object_path(x, step="export_raw_counts", extension="csv")
  if(as_matrix) {
    cnts <- read.csv(cfile, row.names=1)
    return(as.matrix(cnts))
  } else {
    return(read.csv(cfile))
  }
}


#' Get all results of tmod gene set enrichment analysis
#'
#' Get results of the tmod gene set enrichment analysis for all contrasts (default) or a
#' selection of the contrasts.
#' @param x an object of class seasnap_DE_pipeline
#' @param contrasts an optional character vector of the contrasts to return
#' @return a named list; each object of the list is a list of tmod results for the given contrast.
#' @export
get_tmod_res <- function(x, contrasts=NULL) {
  .check_de_obj(x)

  if(is.null(contrasts)) {
    contrasts <- get_contrast_names(x)
  }

  ret <- lapply(contrasts, function(c) {
    sel <- with(x$file_tab, step == "tmod" & contrast == c & extension == "rds")
    path <- x$file_tab$filename[sel]
    seasnap_readRDS(x, path)
  })

  names(ret) <- contrasts
  return(ret)
}

#' Return the tmod gene set database objects
#'
#' Return the tmod gene set database objects
#' @param x an object of class seasnap_DE_pipeline
#' @return a list containing gene set databases and corresponding tmod objects; 
#' @export
get_tmod_dbs <- function(x) {
  get_object(x, step="tmod_dbs", contrast="all", multiple_ok=FALSE)
}

#' Return the tmod database mapping object
#'
#' Return the tmod database mapping object
#'
#' It is best to use the output of this function with [tmod_db_map_ids()].
#'
#' tmod mapping object is a set of mappings between the identifiers used in
#' the tmod databases such as MSigDB and the primary identifiers used in the
#' pipeline. The object is a list consisting of three elements:
#'
#'  * `dbs`: a named character vector. Names are the IDs of the tmod
#'           databases; values are the IDs of the mappings.
#'  * `maps`: a named list of named character vectors. Names of the
#' elements of the list are mapping IDs and correspond to the values of the
#' `dbs` element. Elements of the list are mappings. Each mapping is a
#' named character vector. Names of this vector are PrimaryIDs, and values are
#' the corresponding IDs in the given database.
#' Return the tmod database mapping object
#' @param x an object of class seasnap_DE_pipeline
#' @return a list containing the tmod mapping object
#' @seealso [get_tmod_dbs()], [tmod_db_map_ids()]
#' @export
get_tmod_mapping <- function(x) {
  get_object(x, step="tmod_dbs", extension="mapping.rds", contrast="all", multiple_ok=FALSE)
}


#' Get the DESeq2 object
#'
#' Get the DESeq2 object created by the DE pipeline
#' @param x an object of class seasnap_DE_pipeline
#' @return an object of class DESeq2
#' @export
get_deseq2 <- function(x) {
  get_object(x, step="DESeq2", extension="deseq2.rds", contrast="all", multiple_ok=FALSE)
}

#' Get the covariates
#'
#' Get the covariates from the DESeq2 object
#'
#' Note: rather than reading the covariate file (which might have changed
#' since DESeq2 was run) the function reads the DESeq2 object and gets the
#' column data of that object.
#'
#' Note 2: The returned object is a color data frame. Some columns are 
#' hidden on purpose because they clutter the view. Use `hide=Inf` to
#' prevent this behavior.
#'
#' Note 3: All factors are converted back to strings.
#' @param x an object of class seasnap_DE_pipeline
#' @param hide hide all columns wider than this
#' @return a color data frame with the covariates
#' @export
get_covariates <- function(x, hide=75) {
  .check_de_obj(x)

  ds2 <- get_deseq2(x)
  ret <- as.colorDF(as.data.frame(ds2@colData))
  #ret <- as.data.frame(ds2@colData)

  col_type(ret, c("md5", "filename")) <- "hidden"
  max_w <- sapply(ret, function(x) max(nchar(as.character(x))))
  tohide <- colnames(ret)[max_w >= 75]
  col_type(ret, tohide) <- "hidden"

  isfac <- sapply(ret, function(x) "factor" %in% class(x))

  for(cn in colnames(ret)[isfac]) {
    ret[[cn]] <- as.character(ret[[cn]])
  }

  return(ret)
}

#' Get config
#'
#' Get the configuration of pipeline read from yaml file produced in the pipeline_report step.
#' 
#' @param x an object of class seasnap_DE_pipeline
#' @return a list object.
#' @export
get_config <- function(x) {
  .check_de_obj(x)

  return(x$config)
}



#' Get the pipeline annotation 
#'
#' Get the annotation data frame of the pipeline
#' @param x an object of class seasnap_DE_pipeline
#' @return data frame with the annotation
#' @export
get_annot <- function(x) {
  .check_de_obj(x)
  step      <- "annotation"
  extension <- "rds"

  get_object(x, step, multiple_ok=FALSE)
}

#' Get a pipeline object
#'
#' Get a pipeline object from a given step and extension
#' @param x an object of class seasnap_DE_pipeline
#' @param step of the pipeline
#' @param extension of the object
#' @param multiple_ok Whether it is OK to return multiple objects
#' @param as_list Always return a list of objects
#' @param contrast optional: choice of the contrast
#' @return an object return by readRDS, or a list of objects
#' @export
get_object <- function(x, step, extension="rds", contrast=NULL, multiple_ok=TRUE, as_list=FALSE) {

  paths <- get_object_path(x, step, extension, contrast=contrast, multiple_ok=multiple_ok)

  if(length(paths) > 1 || as_list) {
    ret <- lapply(paths, function(p) seasnap_readRDS(x, p))
  } else {
    ret <- seasnap_readRDS(x, paths)
  }

  return(ret)
}


#' Get file path of a pipeline object
#'
#' Get file path of a pipeline object
#' @param x an object of class seasnap_DE_pipeline
#' @param step of the pipeline
#' @param extension of the object
#' @param contrast optional: choice of the contrast
#' @param multiple_ok whether it is OK to return multiple objects
#' @return a character vector of file paths
#' @importFrom glue glue
#' @export
get_object_path <- function(x, step, extension, contrast=NULL, multiple_ok=TRUE) {

  .check_de_obj(x)
  ft <- x$file_tab
  sel <- ft$step == step
  if(!any(sel)) {
    stop(glue("Step {step} not found"))
  }

  sel <- sel & ft$extension == extension
  if(!any(sel)) {
    stop(glue("Step {step} does not have objects with extension {extension}"))
  }

  if(!is.null(contrast)) {
    sel <- sel & ft$contrast %in% contrast
    if(!any(sel)) {
      contrast <- paste(contrast, ", ")
      stop(glue("No objects for this combination: step={step}, extension={extension}, contrast={contrast}"))
    }
  }

  if(sum(sel, na.rm=TRUE) > 1 & !multiple_ok) {
    stop(glue("More than 1 object in step {step} with extension {extension}"))
  }

  path <- x$file_tab$filename[sel]
  names(path) <- x$file_tab$contrast[sel]
  return(path)
}

#' Show all the steps of the pipeline
#'
#' Show the steps the pipeline has
#' @param x an object of class seasnap_DE_pipeline
#' @return a character vector with the steps
#' @export
get_steps <- function(x) {
  .check_de_obj(x)

  sort(unique(x$file_tab$step))
}



#' Show all the extensions of the pipeline
#'
#' Show the extensions (file types) the pipeline has
#' @param x an object of class seasnap_DE_pipeline
#' @param steps a character vector showing the steps for which the
#'        extensions are to be shown
#' @return a character vector with the extensions
#' @export
get_extensions <- function(x, steps=NULL) {
  .check_de_obj(x)

  if(!is.null(steps)) {
    sel <- x$file_tab$step %in% steps
  } else {
    sel <- 1:nrow(x$file_tab)
  }

  sort(unique(x$file_tab[sel,,drop=FALSE]$extension))
}





