# some packages break in R 4
FROM rocker/tidyverse:4.0.5

LABEL maintainer="tyrone.chen@monash.edu"

RUN R -e 'install_me <- c("argparser", "brio", "colorspace", "diffobj", "dplyr", "ellipsis", "farver", "ggplot2", "ggrepel", "igraph", "isoband", "matrixStats", "mixOmics", "parallel", "plyr", "rARPACK", "Rcpp", "RcppEigen", "reshape2", "rjson", "RSpectra", "stringi", "testthat", "tibble", "tidyr", "utf8", "vctrs", "zeallot"); sapply(install_me, install.packages); devtools::install_github("mixOmicsTeam/mixOmics@devel"); devtools::install_gitlab("tyagilab/sars-cov-2", subdir="multiomics")'
