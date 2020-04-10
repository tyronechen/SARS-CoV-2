#!/usr/bin/R
# combine translatome and proteomics data for sars-cov-2
# data originally from DOI:10.21203/rs.3.rs-17218/v1 - supp tables 1 and 2
library(argparser, quietly=TRUE)
library(ggplot2)
library(parallel)
library(reshape2)
source(file="multiomics_sars-cov-2.R")

parse_argv = function() {
  library(argparser, quietly=TRUE)
  p = arg_parser("Assess the quality of individual omics data modality blocks, \
  Imputes data, Run DIABLO on multi-omics data. Take tsv file of classes, \
  tsv files of omics data (at least 2!) and identify correlation between \
  features. For more information please refer to mixOmics and case studies at \
  http://mixomics.org/mixdiablo/case-study-tcga/")

  # Add command line arguments
  p = add_argument(p, "classes", help="sample information", type="character")
  p = add_argument(p, "--classes_secondary", type="character", default=NA,
    help="secondary sample information (same format as classes)"
  )
  p = add_argument(p, "--data", type="character", nargs=Inf,
    help="paths to omics data. names format: SAMPLEID_OMICTYPE_OPTIONALFIELDS"
  )
  p = add_argument(p, "--dropna_classes", type="character", default=TRUE,
    help="where all replicates for >= 1 class are NA, drop that feature. \
    If both --dropna_classes and --dropna_prop are enabled, perform \
    --dropna_classes and then --dropna_prop."
  )
  p = add_argument(p, "--dropna_prop", type="integer", default=0,
    help="drop feature that does not meet NA proportion threshold, eg if 0.3, \
    drop a feature if >= 0.3 of values are NA. If both --dropna_classes and \
    --dropna_prop are enabled, perform --dropna_classes and then --dropna_prop."
  )
  p = add_argument(p, "--ncpus", help="number of cpus", type="integer", default=2)
  p = add_argument(p, "--dcomp", type="integer", default=0,
    help="number of diablo components (set manually if you get inference error)"
  )
  p = add_argument(p, "--icomp", type="integer", default=10,
    help="component number for imputing (set 0 for no imputation)"
  )
  p = add_argument(p, "--rdata", type="character", default="./diablo.RData",
    help="write RData object here, has (classes, data, diablo, mdist)"
  )
  p = add_argument(p, "--plot", type="character", default="./Rplots.pdf",
    help="write R plots here (will overwrite existing!)"
  )
  p = add_argument(p, "--pcomp", type="integer", default=0,
    help="number of principal components (defaults to number of samples)"
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

  # parse out identifiers coded within the file paths (hardcoded)
  names = sapply(sapply(lapply(paths, strsplit, "/"), tail, 1), tail, 1)
  names = unname(lapply(sapply(names, strsplit, ".", fixed=TRUE), head, 1))
  names = unname(sapply(sapply(names, head, 1), strsplit, "_"))
  names = unlist(lapply(lapply(names, tail, -1), paste, collapse="_"))
  print("Omics data types")
  print(names)

  print(paste("Saving plots to (overwriting existing):", argv$plot))
  pdf(argv$plot)

  # load data and drop features / columns with all NA
  data = lapply(paths, parse_data, missing_as=NA, rmna=TRUE)
  names(data) = names

  # show proportion of NA values in unfiltered data
  mapply(function(x, y) show_na_prop(x, y), data, names)

  # drop features / columns where >= 1 class is not represented
  if (argv$dropna_classes == TRUE) {
    data = lapply(data, remove_na_class, classes)
    save(classes, data, mdist, file=argv$rdata)
  }

  # drop features / columns >= a threshold of NA values
  if (argv$dropna_prop > 0) {
    data = remove_na_prop(data, classes, pch=pch, na_prop=argv$dropna_prop)
    save(classes, data, mdist, file=argv$rdata)
  }

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

  # count missing data after all filtering
  missing = lapply(data, count_missing)
  pca_withna = plot_pca_single(
    data, classes, pch=pch, ncomp=argv$pcomp,
    title=paste("With NA. PC:", argv$pcomp)
  )
  save(classes, data, pca_withna, mdist, file=argv$rdata)

  # impute data if components given
  # refer to http://mixomics.org/methods/missing-values/
  print(argv$icomp)
  if (argv$icomp > 0) {
    print("Impute components set, imputing NA values (set -i 0 to disable)")
    data_imp = impute_missing(data, rep(argv$icomp, length(data)))
    data = replace_missing(data, data_imp)
    pca_impute = plot_pca_single(
      data_imp, classes, pch=pch, ncomp=argv$pcomp,
      title=paste("Imputed. PC:", argv$pcomp, "IC:", argv$icomp)
    )
    save(classes, data, pca_withna, pca_impute, mdist, file=argv$rdata)
    heatmaps = cor_imputed_unimputed(pca_withna, pca_impute, names)
  } else {
    print("Impute components unset, not imputing NA (set -i > 0 to enable)")
    pca_impute = NA
  }

  # multilevel decomposition if secondary variables are specified
  # refer to http://mixomics.org/case-studies/multilevel-vac18/
  if (!is.na(pch)) {
    plot_pca_multilevel(data, classes, pch=pch, ncomp=argv$pcomp,
      title=paste("With NA. PC:", argv$pcomp)
    )
    if (exists("data_imp")) {
      plot_pca_multilevel(data_imp, classes, pch=pch, ncomp=argv$pcomp,
        title=paste("Imputed. PC:", argv$pcomp, "IC:", argv$icomp)
      )
    }
  }

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
  save(classes, data, pca_withna, pca_impute, mdist, file=argv$rdata)

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
  save(classes, data, diablo, pca_withna, pca_impute, mdist, file=argv$rdata)
  print(paste("Saved plots to:", argv$plot))
  dev.off()
}

main()
