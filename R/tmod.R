#' Map the PrimaryIDs for a database
#'
#' Map the PrimaryIDs for a database
#'
#' Using precomputed mappings, map PrimaryIDs (usually ENSEMBL IDs) to the
#' IDs of the desired database.
#' @param ids character vector of PrimaryIDs
#' @param tmod_dbs_mapping_obj the object returned by `get_object(x, "tmod_dbs", "mapping.rds")`
#' @return a character vector
#' @export
tmod_db_map_ids <- function(x, ids, dbname, tmod_dbs_mapping_obj=NULL) {

  .check_de_obj(x)

  if(is.null(tmod_dbs_mapping_obj)) {
    tmod_dbs_mapping_obj <- get_object(x, step="tmod_dbs", extension="mapping.rds")
  }

  if(!dbname %in% names(tmod_dbs_mapping_obj$dbs)) {
    stop(glue("No mapping for db {dbname} in the mapping object"))
  }

  mapping_id <- tmod_dbs_mapping_obj$dbs[dbname]
  mapping    <- tmod_dbs_mapping_obj$maps[[mapping_id]]

  ret <- mapping[ids]
  if(all(is.na(ret))) {
    warning("No IDs were found in the mapping... are you sure you are using PrimaryIDs?\nCheck the annotation data frame (see `get_annot()`)")
  }
  names(ret) <- ids
  return(ret)
}





#' Create tmod evidence plot
#'
#' Create an evidence plot (ROC curve) for the given module and database
#'
#' Evidence plots are ROC curves which visualize the strength of a gene set
#' enrichment based on an ordered list of genes.
#' For a given contrast, given gene set ID belonging to a given database,
#' `plot_evidence` produces such a plot.
#'
#' The contrast parameter can either be a single character value (name of
#' the contrast) or a character vector of gene IDs. In
#' the first case, the argument is the name of the contrast to be used. 
#' In the second case, the vector is assumed to be a vector of PrimaryIDs
#' (typically ENSEMBL IDs), just like the column "PrimaryID" in the
#' annotation data frame (returned by get_annot).
#'
#' Since loading of objects can take time, the larger objects (tmod dbs,
#' contrasts, precomputed gene lists) can be provided as optional named arguments.
#' @param x an object of class seasnap_DE_pipeline
#' @param id ID of the gene set (module) 
#' @param dbname Name of the tmod database (see get_tmod_db_names)
#' @param contrast either a character vector of length one or a character vector of gene IDs
#' @param contrast_obj a data frame with columns named PrimaryID, log2FoldChange and p.adj
#' @param annot_obj annotation data frame returned by get_annot(x)
#' @param tmod_dbs_obj the object returned by `get_object(x, "tmod_dbs")`
#' @param tmod_dbs_mapping_obj the object returned by `get_object(x, "tmod_dbs", "mapping.rds")`
#' @param gl (optional) gene list object returned by get_tmod_gl (if contrast is a character vector of length 1)
#' @param sort which sorting type (must be present in the gene list, don't worry if you don't know what that is)
#' @param ... further arguments passed to the tmod::evidencePlot() function
#' @return NULL
#' @import tmod
#' @export
plot_evidence <- function(x, id, dbname, contrast, sort="pval", contrast_obj=NULL, gl=NULL, 
  tmod_dbs_obj=NULL,
  tmod_dbs_mapping_obj=NULL, annot_obj=NULL, ...) {
  .check_de_obj(x)

  if(is.null(annot_obj)) {
    annot_obj <- get_annot(x)
  }

  if(length(contrast) == 1) {
    cn <- get_contrast_names(x)
    if(!contrast %in% cn) {
      cn <- paste(cn, collapse=", ")
      stop(glue("Contrast must be one of {cn}"))
    }

    if(is.null(contrast_obj)) {
      contrast_obj <- get_contrasts(x, contrast)[[1]]
    }

    if(is.null(gl)) {
      gl <- get_object(x, step="tmod", extension="gl.rds", contrast=contrast, multiple_ok=FALSE)
    }
    gl <- gl[[dbname]][[sort]]
    primary_ids <- names(gl)
  } else {
    ## we need to translate the PrimaryIDs
    gl <- tmod_db_map_ids(x, contrast, dbname, tmod_dbs_mapping_obj)
    primary_ids <- contrast
  }

  ## by now, we have the gl and primary_ids defined and mapped to the correct db

  ## get symbols instead of primary ids if possible
  symbols <- primary_ids
  if("SYMBOL" in colnames(annot_obj)) {
    symbols <- annot_obj$SYMBOL[ match(primary_ids, annot_obj$PrimaryID) ]
  }

  ## now for some color
  if(is.null(contrast_obj)) {
    colors <- NULL
  } else {
    if(!all(c("log2FoldChange", "p.adj") %in% colnames(contrast_obj))) {
      stop("Contrast object must have columns log2FoldChange and p.adj")
    }

    if(!any(primary_ids %in% rownames(contrast_obj))) {
      stop("none of the primary IDs is in rownames(contrast_obj)")
    }

    contrast_obj <- contrast_obj[ primary_ids, , drop=FALSE ]
    p <- contrast_obj$p.adj
    l <- contrast_obj$log2FoldChange
    colors <- ifelse(l < 0, 
      ifelse(p < .05, 'blue', '#660000'),
      ifelse(p < .05, 'red', '#000066'))
    names(colors) <- gl
  }

  if(is.null(tmod_dbs_obj)) {
    tmod_dbs_obj <- get_object(x, "tmod_dbs")
  }

  mset <- tmod_dbs_obj[[dbname]][["dbobj"]]





}
