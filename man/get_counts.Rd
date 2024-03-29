% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_utils.R
\name{get_counts}
\alias{get_counts}
\title{Get all counts from the pipeline}
\usage{
get_counts(x, as_matrix = TRUE)
}
\arguments{
\item{x}{an object of class seasnap_DE_pipeline}

\item{as_matrix}{if FALSE, a data frame will be returned}
}
\value{
either a matrix with raw counts, rownames corresponding to PrimaryIDs
and column names corresponding to sample names, or data frame
with first column being the PrimaryIDs
}
\description{
Get counts from the import_gene_counts step
}
