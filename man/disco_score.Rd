% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{disco_score}
\alias{disco_score}
\title{Calculate the disco score}
\usage{
disco_score(
  contrast1,
  contrast2,
  minp = 1e-16,
  by = 0,
  primary_id = "PrimaryID"
)
}
\arguments{
\item{contrast1, contrast2}{data frames with rownames corresponding to
IDs (they don't need to be in the same order) and columns \code{log2FoldChange}
and \code{pvalue}.}

\item{minp}{minimum p-value}

\item{by}{column by which the contrast data frames should be merged
(passed to \code{merge}). Default: merge by row names}

\item{primary_id}{the name which should be assigned to the identifier
column which results from the merge}
}
\value{
a merged data frame containing column "disco.score"
}
\description{
Merge two contrasts and calculate the disco score
}
\details{
Disco score is a heuristic score
reflecting the strength of similarityor dissimilarity between two
log~2~ FC values weighted by the corresponding p-values and given by the
formula log2FC.x * log2FC.y * (-log10(pval.x) - log10(pval.y)).
}
