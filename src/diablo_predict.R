#!/usr/bin/R
# combine translatome and proteomics data for sars-cov-2
# data originally from DOI:10.21203/rs.3.rs-17218/v1 - supp tables 1 and 2
library(argparser, quietly=TRUE)
# library(mixOmics)
source(file="multiomics_sars-cov-2.R")

parse_argv = function() {
  p = arg_parser("Use predicted diablo model")

  # Add command line arguments
  p = add_argument(p, "data", help="load RData object", type="character")

  # Parse the command line arguments
  argv = parse_args(p)

  # Do work based on the passed arguments
  return(argv)
}

main = function() {
  argv = parse_argv()
  load(argv$data)
  plot_diablo(diablo)
  assess_performance(diablo, dist=dist)
  predict_diablo(diablo, data, classes)
}

main()
