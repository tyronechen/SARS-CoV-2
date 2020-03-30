#!/usr/bin/R
# combine translatome and proteomics data for sars-cov-2
# data originally from DOI:10.21203/rs.3.rs-17218/v1 - supp tables 1 and 2
library(argparser, quietly=TRUE)
library(parallel)
source(file="multiomics_sars-cov-2.R")

parse_argv = function() {
  library(argparser, quietly=TRUE)
  p = arg_parser("Run DIABLO on multi-omics data")

  # Add command line arguments
  p = add_argument(p, "classes", help="sample information", type="character")
  p = add_argument(p, "--data", help="paths to omics data", type="character", nargs=Inf)
  p = add_argument(p, "--cpus", help="number of cpus", type="int", default=2)
  p = add_argument(p, "--ncomp", help="number of components", type="int", default=0)
  p = add_argument(p, "--out", help="write RData object here", type="character")
  p = add_argument(p, "--distance", help="distance metric to use [max.dist, centroids.dist, mahalanobis.dist]", type="character", default="max.dist")

  # Parse the command line arguments
  argv = parse_args(p)

  # Do work based on the passed arguments
  return(argv)
}

main = function() {
  argv = parse_argv()

  print("Available cpus:")
  print(detectCores())
  print("Using cpus (change with --cpus):")
  print(argv$cpus)
  # q()
  print("Distance measure:")
  print(argv$distance)

  options(warn=1)

  paths = argv$data

  print("Paths to data:")
  print(paths)

  print("Parsing classes")
  classes = parse_classes(argv$classes)

  # parse out identifiers coded within the file paths
  names = sapply(sapply(lapply(paths, strsplit, "/"), tail, 1), tail, 1)
  names = unname(lapply(sapply(names, strsplit, ".", fixed=TRUE), head, 1))
  names = unname(sapply(sapply(names, head, 1), strsplit, "_"))
  names = unlist(lapply(lapply(names, tail, -1), paste, collapse="_"))
  print("Omics data types (names follow SAMPLEID_OMICTYPE_OPTIONALFIELDS):")
  print(names)

  data = lapply(paths, parse_data)
  names(data) = names
  print("Data dimensions:")
  dimensions = lapply(data, dim)
  print(dimensions)

  design = create_design(data)

  # check dimension
  print(summary(classes))
  print("Y (classes):")
  print(classes)
  print("Design:")
  print(design)

  plot_individual_blocks(data, classes)

  # NOTE: if you get tuning errors, set ncomp manually with --ncomp N
  if (argv$ncomp == 0) {
    tuned = tune_ncomp(data, classes, design)
    print("Parameters with lowest error rate:")
    tuned = tuned$choice.ncomp$WeightedVote["Overall.BER",]
    ncomp = tuned[which.max(tuned)]
  } else {
    ncomp = argv$ncomp
  }
  print("Number of components:")
  print(ncomp)

  keepx = tune_keepx(data, classes, ncomp, design, cpus=argv$cpus, dist=distance)
  print("keepx:")
  print(keepx)
  diablo = run_diablo(data, classes, ncomp, keepx, design)
  print("diablo design:")
  print(diablo$design)
  # selectVar(diablo, block = "proteome", comp = 1)$proteome$name
  plot_diablo(diablo)
  assess_performance(diablo, dist=distance)
  predict_diablo(diablo, data, classes)

  if (! exists(argv$out)) {
    print(paste("Saving diablo data to:", argv$out))
    save(classes, data, diablo, distance, file=argv$out)
  } else {
    print(paste("Saving diablo data to:", "./diablo.RData"))
    save(classes, data, diablo, distance, file="./diablo.RData")
  }
}

main()
