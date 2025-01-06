
# Rseasnap

<!-- badges: start -->
<!-- badges: end -->

A collection of functions which make it easy to access the files and
objects created by the sea-snap DE pipeline.

Why Rseasnap:

 * shortcuts for accessing pipeline objects
 * produce plots on the fly
 * easier to work with multiple pipelines 

Why not part of sea-snap distribution (e.g. as R scripts):

 * better documentation
 * can be used for other tasks as well
 * easier to maintain

## Installation

You can install the current version of Rseasnap as follows:

``` r
library(devtools)
install_github("bihealth/Rseasnap")
```

Alternatively, clone the git repo and install from the clone directory,
e.g.

``` r
system("git clone git@github.com:bihealth/Rseasnap.git")
devtools::install_git("Rseasnap")
```

## Quickstart

Load the DE pipeline object:

``` r
library(Rseasnap)
pip <- load_de_pipeline("DE_config.yaml")
```

Once the YAML file is loaded into the object, you can use it to access
various objects and files:

``` r
# get the counts
counts <- get_counts(pip)

# get the expression data
expr <- get_exprs(pip)

# show contrast names
get_contrast_names(pip)

# get the DE results
res <- get_contrasts(pip)

# get the gene annotation file
annot <- get_annot(pip)

# get the covariate file
covar <- get_covariates(pip)
```










