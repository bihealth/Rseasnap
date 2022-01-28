
#' @param pc the ID of the principal component to show the evidence plot for
#' @param direction abs, up or down. "abs" means the component weights are
#' sorted descendingly by absolute value; "up" means they are sorted with
#' the smallest at the top of the list; "down" means they are sorted with
#' the largest values at the top of the list.
#' @rdname plot_evidence
#' @export
plot_tmod_pca_evidence <- function (x, id, dbname, pc, direction="abs", tmod_dbs_obj = NULL, 
  annot_obj = NULL, ...) {
  browser()

  .check_de_obj(x)
  if (is.null(tmod_dbs_obj)) {
    message("Loading tmod dbs object (to speed up, use the tmod_dbs_obj option)")
    tmod_dbs_obj <- get_tmod_dbs(x)
  }
  if (is.null(annot_obj)) {
    annot_obj <- get_annot(x)
  }

  if (length(id) != 1L) {
    stop("id must be of length 1")
  }

  if(length(pc) != 1L) {
    stop("only 1 PC can be used")
  }

  tppp <- get_tmod_pca_prcomp(x)
  cn <- colnames(tppp$rotation)

  if (!pc %in% cn) {
    cn <- paste(cn, collapse = ", ")
    stop(glue("PC must be one of {cn}"))
  }

  if (direction == "up"){
    gl <- tppp$rotation[order(tppp$rotation[,pc]),]
  }

  if (direction == "down"){
    gl <- tppp$rotation[order(-tppp$rotation[,pc]),]
  }

  if (direction == "abs"){
    gl <- tppp$rotation[order(-abs(tppp$rotation[,pc])),]
  }

  primary_ids <- rownames(gl)
  gl <- tmod_db_map_ids(x, primary_ids, dbname)

  symbols <- primary_ids
  if ("SYMBOL" %in% colnames(annot_obj)) {
    symbols <- annot_obj$SYMBOL[match(primary_ids, annot_obj$PrimaryID)]
  }
  names(symbols) <- primary_ids
  mset <- tmod_dbs_obj[[dbname]][["dbobj"]]
  genes <- getModuleMembers(id, mset)[[id]]
  if (is.null(genes)) {
    stop(glue("Gene set {id} not found in db {dbname}"))
  }
  sel <- gl %in% genes
  gene.labels <- symbols[sel]
  names(gene.labels) <- gl[sel]
  if (is.null(tmod_dbs_obj)) {
    tmod_dbs_obj <- get_object(x, "tmod_dbs")
    tmod_dbs_obj <- tmod_dbs_obj[dbname]
    if (is.null(tmod_dbs_obj[[dbname]])) {
      stop(glue("No database with name {dbname}"))
    }
  }
  evidencePlot(gl, m = id, mset = mset, 
    gene.labels = gene.labels, ...)
  ret <- list( annot_obj = annot_obj, 
    gl = gl, tmod_dbs_obj = tmod_dbs_obj, 
    x = x, dbname = dbname)
  class(ret) <- c(class(ret), "evidence_plot_data")
  return(invisible(ret))
}
