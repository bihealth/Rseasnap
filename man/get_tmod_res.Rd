% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_utils.R
\name{get_tmod_res}
\alias{get_tmod_res}
\title{Get all results of tmod gene set enrichment analysis}
\usage{
get_tmod_res(x, contrasts = NULL)
}
\arguments{
\item{x}{an object of class seasnap_DE_pipeline}

\item{contrasts}{an optional character vector of the contrasts to return}
}
\value{
a named list; each object of the list is a list of tmod results for the given contrast.
}
\description{
Get results of the tmod gene set enrichment analysis for all contrasts (default) or a
selection of the contrasts.
}
