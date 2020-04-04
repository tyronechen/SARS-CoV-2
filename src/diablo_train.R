#!/usr/bin/R
# combine translatome and proteomics data for sars-cov-2
# data originally from DOI:10.21203/rs.3.rs-17218/v1 - supp tables 1 and 2
library(argparser, quietly=TRUE)
library(parallel)
source(file="multiomics_sars-cov-2.R")

parse_argv = function() {
  library(argparser, quietly=TRUE)
  p = arg_parser("Run DIABLO on multi-omics data. Take tsv file of classes, \
  tsv files of omics data (at least 2!) and identify correlation between \
  features. For more information refer to \
  http://mixomics.org/mixdiablo/case-study-tcga/")

  # Add command line arguments
  p = add_argument(p, "classes", help="sample information", type="character")
  p = add_argument(p, "--classes_secondary", type="character", default=NA,
    help="secondary sample information (same format as classes)"
  )
  p = add_argument(p, "--classes_dropna", type="character", default=TRUE,
    help="where all replicates for >= 1 class are unreported, drop that feature"
  )
  p = add_argument(p, "--data", type="character", nargs=Inf,
    help="paths to omics data. names format: SAMPLEID_OMICTYPE_OPTIONALFIELDS"
  )
  p = add_argument(p, "--ncpus", help="number of cpus", type="int", default=2)
  p = add_argument(p, "--dcomp", type="int", default=0,
    help="number of diablo components (set manually if you get inference error)"
  )
  p = add_argument(p, "--icomp", type="int", default=10,
    help="component number for imputing (set 0 for no imputation)"
  )
  p = add_argument(p, "--rdata", type="character", default="./diablo.RData",
    help="write RData object here, has (classes, data, diablo, mdist)"
  )
  p = add_argument(p, "--plot", help="write R plots here", type="character",
    default="./Rplots.pdf"
  )
  p = add_argument(p, "--mdist", type="character", default="max.dist",
    help="distance metric to use [max.dist, centroids.dist, mahalanobis.dist]"
  )

  # Parse the command line arguments
  argv = parse_args(p)

  # Do work based on the passed arguments
  return(argv)
}

main = function() {
  argv = parse_argv()

  # print some diagnostics for debugging
  print("Available cpus:")
  print(detectCores())
  print("Using cpus (change with --ncpus):")
  print(argv$ncpus)
  print("Distance measure:")
  print(argv$mdist)
  mdist = argv$mdist

  options(warn=1)

  paths = argv$data
  print("Paths to data:")
  print(paths)
  print("Parsing classes")
  classes = parse_classes(argv$classes)

  if (!is.na(argv$classes_secondary)) {
    print("Parsing secondary classes")
    pch = parse_classes(argv$classes_secondary)
  } else {
    pch = NA
  }

  # parse out identifiers coded within the file paths
  names = sapply(sapply(lapply(paths, strsplit, "/"), tail, 1), tail, 1)
  names = unname(lapply(sapply(names, strsplit, ".", fixed=TRUE), head, 1))
  names = unname(sapply(sapply(names, head, 1), strsplit, "_"))
  names = unlist(lapply(lapply(names, tail, -1), paste, collapse="_"))
  print("Omics data types (names follow SAMPLEID_OMICTYPE_OPTIONALFIELDS):")
  print(names)

  # load data and drop cols with all NA
  data = lapply(paths, parse_data, missing_as=NA, rmna=TRUE)
  names(data) = names

  # drop cols where >= 1 class is not represented
  if (argv$classes_dropna == TRUE) {
    data = lapply(data, remove_na_class, classes)
  }
  save(classes, data, mdist, file=argv$rdata)

  # check dimensions
  print("Data dimensions:")
  dimensions = lapply(data, dim)
  print(dimensions)

  design = create_design(data)

  # check classes
  print(summary(classes))
  print("Y (classes):")
  print(classes)
  print("Design:")
  print(design)
  print(paste("Saving plots to:", argv$plot))
  pdf(argv$plot)

  # count missing data
  missing = lapply(data, count_missing)
  data_pca = plot_individual_blocks(data, classes, pch=pch, title="Not imputed")

  # impute data if components given
  print(argv$icomp)
  if (argv$icomp > 0) {
    print("Impute components set, imputing NA values (set -i 0 to disable)")
    data_imp = impute_missing(data, rep(argv$icomp, length(data)))
    data = replace_missing(data, data_imp)
    data_pca = plot_individual_blocks(data, classes, pch=pch, title="Imputed")
  }

  # plot pcas for each block
  # save(classes, data, mdist, file=argv$rdata)
  # data_pca = plot_individual_blocks(data, classes, pch, title="Imputed")
  save(classes, data, data_pca, mdist, file=argv$rdata)

  # NOTE: if you get tuning errors, set dcomp manually with --dcomp N
  if (argv$dcomp == 0) {
    tuned = tune_ncomp(data, classes, design)
    print("Parameters with lowest error rate:")
    tuned = tuned$choice.ncomp$WeightedVote["Overall.BER",]
    dcomp = tuned[which.max(tuned)]
  } else {
    dcomp = argv$dcomp
  }
  print("Number of components:")
  print(dcomp)

  # remove invariant columns
  data = lapply(data, remove_novar)
  save(classes, data, data_pca, mdist, file=argv$rdata)

  # tune diablo parameters and run diablo
  keepx = tune_keepx(data, classes, dcomp, design, cpus=argv$ncpus, dist=mdist)
  print("keepx:")
  print(keepx)
  diablo = run_diablo(data, classes, dcomp, keepx, design)
  print("diablo design:")
  print(diablo$design)
  # selectVar(diablo, block = "proteome", comp = 1)$proteome$name
  plot_diablo(diablo)
  assess_performance(diablo, dist=mdist)
  predict_diablo(diablo, data, classes)

  # save RData object for future reference
  print(paste("Saving diablo data to:", argv$rdata))
  save(classes, data, data_pca, diablo, mdist, file=argv$rdata)
  print(paste("Saved plots to:", argv$plot))
  dev.off()
}

main()
