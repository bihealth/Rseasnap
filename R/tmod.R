#' Map the PrimaryIDs for a database
#'
#' Map the PrimaryIDs for a database
#'
#' Using precomputed mappings, map PrimaryIDs (usually ENSEMBL IDs) to the
#' IDs of the desired database.
#' @param x an object of class seasnap_DE_pipeline
#' @param ids character vector of PrimaryIDs
#' @param tmod_dbs_mapping_obj the object returned by `get_object(x, "tmod_dbs", "mapping.rds")`
#' @param dbname string, name of the database to use
#' @return a named character vector of the same length and order as `ids`
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
#' It is usually much more efficient to get the tmod database object first
#' and pass it as an optional argument than to load it each time this function
#' is called. The mapping object (tmod_dbs_mapping_obj) is less critical,
#' as it is usually much smaller and faster to load.
#'
#' Alternatively, if you assign the output of `plot_evidence()` to a
#' variable, this variable will store the necessary object and can be used
#' with the `plot()` function to efficiently plot the data again.
#' @param x an object of class seasnap_DE_pipeline
#' @param id ID of the gene set (module) 
#' @param dbname Name of the tmod database (see get_tmod_db_names)
#' @param contrast either a character vector of length one or a character vector of gene IDs
#' @param contrast_obj a data frame with columns named PrimaryID, log2FoldChange and padj
#' @param annot_obj annotation data frame returned by get_annot(x)
#' @param tmod_dbs_obj the object returned by `get_object(x, "tmod_dbs")`
#' @param tmod_dbs_mapping_obj the object returned by `get_object(x, "tmod_dbs", "mapping.rds")`
#' @param gl (optional) gene list object returned by get_tmod_gl (if contrast is a character vector of length 1)
#' @param sort which sorting type (must be present in the gene list, don't worry if you don't know what that is)
#' @param ... further arguments passed to the tmod::evidencePlot() function
#' @return Returns invisibly an object of class `evidence_plot_data` which
#' can be reused to quickly repeat the plot with another gene set ID using the
#' `plot.evidence_plot_data` function.
#' @import tmod
#' @export
plot_evidence <- function(x, id, dbname, contrast, sort="pval", contrast_obj=NULL, gl=NULL, 
  tmod_dbs_obj=NULL,
  tmod_dbs_mapping_obj=NULL, annot_obj=NULL, ...) {

  .check_de_obj(x)

  if(is.null(annot_obj)) { annot_obj <- get_annot(x) }
  if(length(id) != 1L) {
    stop("id must be of length 1")
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
      gl <- gl[[dbname]][[sort]]
    }
    primary_ids <- names(gl)
  } else {
    ## we need to translate the PrimaryIDs
    gl <- tmod_db_map_ids(x, contrast, dbname, tmod_dbs_mapping_obj)
    primary_ids <- contrast
  }

  ## by now, we have the gl and primary_ids defined and mapped to the correct db

  ## get symbols instead of primary ids if possible
  symbols <- primary_ids
  if("SYMBOL" %in% colnames(annot_obj)) {
    symbols <- annot_obj$SYMBOL[ match(primary_ids, annot_obj$PrimaryID) ]
  }
  names(symbols) <- primary_ids

  ## We need to figure out what genes are in the gene set to select them as
  ## labels
  mset  <- tmod_dbs_obj[[dbname]][["dbobj"]]
  genes <- getModuleMembers(id, mset)[[id]]
  if(is.null(genes)) {
    stop(glue("Gene set {id} not found in db {dbname}"))
  }

  sel                <- gl %in% genes
  gene.labels        <- symbols[sel]
  names(gene.labels) <- gl[sel]

  ## now for some color
  if(is.null(contrast_obj)) {
    colors <- NULL
  } else {
    if(!all(c("log2FoldChange", "padj") %in% colnames(contrast_obj))) {
      stop("Contrast object must have columns log2FoldChange and padj")
    }

    if(!any(primary_ids %in% rownames(contrast_obj))) {
      stop("none of the primary IDs is in rownames(contrast_obj)")
    }

    contrast_obj <- contrast_obj[ primary_ids, , drop=FALSE ]
    p <- contrast_obj$padj
    l <- contrast_obj$log2FoldChange
    colors <- ifelse(l < 0, 
      ifelse(p < .05, 'blue', '#660000'),
      ifelse(p < .05, 'red', '#000066'))
    names(colors) <- gl
  }

  if(is.null(tmod_dbs_obj)) {
    tmod_dbs_obj <- get_object(x, "tmod_dbs")
    tmod_dbs_obj <- tmod_dbs_obj[dbname]
    if(is.null(tmod_dbs_obj[[dbname]])) {
      stop(glue("No database with name {dbname}"))
    }
  }

  evidencePlot(gl, m=id, mset=mset, gene.colors=colors, gene.labels=gene.labels, ...)

  ret <- list(contrast_obj=contrast_obj, 
              annot_obj=annot_obj,
              gl=gl,
              contrast=contrast,
              tmod_dbs_obj=tmod_dbs_obj,
              x=x,
              dbname=dbname)
  class(ret) <- c(class(ret), "evidence_plot_data")
  return(invisible(ret))

}

#' @rdname plot_evidence
#' @export
plot.evidence_plot_data <- function(x, id, ...) {

  if(!"evidence_plot_data" %in% class(x)) {
    stop("x must be of class evidence_plot_data")
  }

  args    <- x
  args$id <- id
  do.call(plot_evidence, args)
}

#' @export
print.evidence_plot_data <- function(x, ...) {

  cat(glue(
"Object of class evidence_plot_data for gene set database {x$dbname}.
Use plot(object, id) to show evidence plot
for other gene set IDs for the same database and contrast\n\n"))
}




#' Run tmod gene set enrichment
#'
#' Run tmod gene set enrichment
#'
#' It is usually much more efficient to get the tmod database object first
#' and pass it as an optional argument than to load it each time this function
#' is called. The mapping object (tmod_dbs_mapping_obj) is less critical,
#' as it is usually much smaller and faster to load.
#' @param x an object of class seasnap_DE_pipeline
#' @param gl an ordered character vector with primary IDs
#' @param dbname name of the gene set database to use
#' @param tmod_dbs_obj object return by the `get_tmod_dbs()` function
#' @param tmod_dbs_mapping_obj the object returned by `get_object(x, "tmod_dbs", "mapping.rds")`
#' @param func tmod gene set enrichment function (by default, tmodCERNOtest)
#' @param ... further arguments passed to the tmod gene set enrichment test function
#' @import tmod
#' @seealso [tmod::tmodCERNOtest()], [get_tmod_dbs()], [get_tmod_mapping()]
#' @export
test_gsea_tmod <- function(x, gl, dbname, tmod_dbs_obj=NULL, tmod_dbs_mapping_obj=NULL, func=tmodCERNOtest, ...) {

  if(is.null(tmod_dbs_mapping_obj)) {
    tmod_dbs_mapping_obj <- get_object(x, step="tmod_dbs", extension="mapping.rds")
    if(is.na(tmod_dbs_mapping_obj$dbs[dbname])) {
      stop(glue("No mapping for database {dbname} in the mapping object"))
    }
  }

  if(is.null(tmod_dbs_obj)) {
    tmod_dbs_obj <- get_tmod_dbs(x)[dbname]
    if(is.null(tmod_dbs_obj[[dbname]])) {
      stop(glue("No database with name {dbname}"))
    }
  }

  gl   <- tmod_db_map_ids(x, gl, dbname, tmod_dbs_mapping_obj)
  mset <- tmod_dbs_obj[[dbname]][["dbobj"]]

  args <- c(list(l=gl, mset=mset), list(...))
  do.call(func, args)
}

#' Produce a tabbed DT results table
#'
#' Produce a tabbed DT results table of gene set enrichment test.
#' @param res results of a gene set enrichment test (e.g.  [test_gsea_tmod()])
#' @param pval_thr report only results with FDR < pval_thr
#' @param es_trh report only results with effect size > es_thr
#' @import DT
#' @importFrom knitr knit_child
#' @seealso [test_gsea_tmod()]
#' @export
test_results_tabbed_DT <- function(res, pval_thr=.05, es_thr=.6) {

  ## table template
  template <- "
\n```{r}
datatable(.r, extensions=c('Buttons','FixedColumns'), rownames=FALSE, escape=FALSE,
        options = list(scrollX=TRUE, fixedColumns=list(leftColumns=1), dom='Bfrtip', buttons=c('excel', 'csv'))) %>%
        formatSignif('p.value', 2) %>%
        formatSignif('FDR', 2) %>%
        formatSignif('AUC', 2)

\n```

"

  cat("\n##### {.tabset}\n")

  for(db in names(res)) {
    cat(sprintf("\n###### %s\n", db))
    .r <- res[[db]]
    if(is.null(.r)) {
      cat("\nNo results.\n")
    } else {
      .r <- .r %>% dplyr::filter(AUC > es_thr & adj.P.Val < pval_thr)
      if(nrow(.r) < 1) {
        cat(sprintf("\nNo results below specified thresholds (ES > %.2f, p_val < %.2f).\n", es_thr, pval_thr))
      } else {
        out <- knitr::knit_child(text=template, quiet=TRUE)
        cat(out)
      }
    }
  }
  cat("\n##### {-}\n")
}



