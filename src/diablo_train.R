#!/usr/bin/R
# combine translatome and proteomics data for sars-cov-2
# data originally from DOI:10.21203/rs.3.rs-17218/v1 - supp tables 1 and 2
library(argparser, quietly=TRUE)
# library(mixOmics)
source(file="multiomics_sars-cov-2.R")

parse_argv = function() {
  library(argparser, quietly=TRUE)
  p = arg_parser("Run DIABLO on multi-omics data")

  # Add command line arguments
  p = add_argument(p, "proteome", help="proteome path", type="character")
  p = add_argument(p, "translatome", help="translatome path", type="character")
  p = add_argument(p, "classes", help="classes path", type="character")
  p = add_argument(p, "--cpus", help="number of cpus", type="int", default=2)
  p = add_argument(p, "--out", help="write RData object", type="character")
  p = add_argument(p, "--dist", help="distance metric to use", type="character")

  # Parse the command line arguments
  argv = parse_args(p)

  # Do work based on the passed arguments
  return(argv)
}

main = function() {
  argv = parse_argv()
  if (! exists(argv$dist)) {
    dist = argv$dist
  } else {
    dist = "centroids.dist"
  }
  print("Distance measure:")
  print(dist)
  options(warn=1)
  # "../data/proteome_diablo.txt"
  # "../data/translatome_diablo.txt"
  # "../data/classes_diablo.txt"
  prot = parse_data(argv$proteome)# + 1
  tran = parse_data(argv$translatome)# + 1
  print(dim(prot))
  print(dim(tran))
  classes = parse_classes(argv$classes)
  data = list(proteome = prot, translatome = tran)
  design = create_design(data)

  # check dimension
  print(summary(classes))
  print("Y (classes):")
  print(classes)
  print("Design:")
  print(design)

  # NOTE: if you get tuning errors, disable this block and set ncomp manually
  tuned = tune_ncomp(data, classes, design)
  print("Parameters with lowest error rate:")
  tuned = tuned$choice.ncomp$WeightedVote["Overall.BER",]
  ncomp = tuned[which.max(tuned)]

  # ncomp = length(unique(classes))
  # ncomp = 10
  print("Components:")
  print(ncomp)
  keepx = tune_keepx(data, classes, ncomp, design, cpus=argv$cpus, dist=dist)
  print("keepx:")
  print(keepx)
  diablo = run_diablo(data, classes, ncomp, keepx, design)
  print("diablo design:")
  print(diablo$design)
  # selectVar(diablo, block = "proteome", comp = 1)$proteome$name
  plot_diablo(diablo)
  assess_performance(diablo, dist=dist)
  predict_diablo(diablo, data, classes)

  if (! exists(argv$out)) {
    print(paste("Saving diablo data to:", argv$out))
    save(classes, data, diablo, dist, file=argv$out)
  } else {
    print(paste("Saving diablo data to:", "./diablo.RData"))
    save(classes, data, diablo, dist, file="./diablo.RData")
  }
}

main()
