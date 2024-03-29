% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmod.R
\name{tmod_db_map_ids}
\alias{tmod_db_map_ids}
\alias{tmod_rev_db_map_ids}
\title{Map the PrimaryIDs for a database}
\usage{
tmod_db_map_ids(x, ids, dbname, tmod_dbs_mapping_obj = NULL)

tmod_rev_db_map_ids(x, ids, dbname, tmod_dbs_mapping_obj = NULL)
}
\arguments{
\item{x}{an object of class seasnap_DE_pipeline}

\item{ids}{character vector of PrimaryIDs (for tmod_db_map_ids) or tmod db ids (for tmod_rev_db_map_ids)}

\item{dbname}{string, name of the database to use}

\item{tmod_dbs_mapping_obj}{the object returned by \code{get_object(x, "tmod_dbs", "mapping.rds")}}
}
\value{
a named character vector of the same length and order as \code{ids}
}
\description{
Map the PrimaryIDs for a database
}
\details{
Using precomputed mappings, map PrimaryIDs (usually ENSEMBL IDs) to the
IDs of the desired database (tmod_db_map_ids) or map the IDs from a data
base to PrimaryIDs (usually ENSEMBL IDs; tmod_rev_db_map_ids)
}
