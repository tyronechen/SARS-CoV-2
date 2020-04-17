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
    help="secondary sample information eg individual (same format as classes)"
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
  p = add_argument(p, "--plsdacomp", type="integer", default=0,
    help="number of components for plsda (defaults to number of samples)"
  )
  p = add_argument(p, "--splsdacomp", type="integer", default=0,
    help="number of components for splsda (defaults to number of samples)"
  )
  p = add_argument(p, "--splsda_keepx", type="vector", default=NA, nargs="+",
    help="variables to keep for splsda"
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

  # initialise plots
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
    save(classes, data, mdist, argv, file=argv$rdata)
  }

  # drop features / columns >= a threshold of NA values
  if (argv$dropna_prop > 0) {
    data = remove_na_prop(data, classes, pch=pch, na_prop=argv$dropna_prop)
    save(classes, data, mdist, argv, file=argv$rdata)
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
  save(classes, data, pca_withna, mdist, argv, file=argv$rdata)

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
    heatmaps = cor_imputed_unimputed(pca_withna, pca_impute, names)
  } else {
    print("Impute components unset, not imputing NA (set -i > 0 to enable)")
    data_imp = NA
    pca_impute = NA
  }
  save(classes, data, data_imp, pca_withna, pca_impute, mdist, argv,
    file=argv$rdata
  )

  # multilevel decomposition if secondary variables are specified
  # refer to http://mixomics.org/case-studies/multilevel-vac18/
  # multilevel pca
  if (!is.na(data_imp)) {
    input_data = data_imp
  } else {
    input_data = data
  }

  if (!is.na(pch)) {
    data_pca_multilevel = plot_pca_multilevel(
      input_data, classes, pch=pch, ncomp=argv$pcomp,
      title=paste("With NA. PC:", argv$pcomp)
    )
    if (exists("data_imp")) {
      data_pca_multilevel = plot_pca_multilevel(
        input_data, classes, pch=pch, ncomp=argv$pcomp,
        title=paste("Imputed. PC:", argv$pcomp, "IC:", argv$icomp)
      )
    }
  } else { data_pca_multilevel = NA }
  save(classes, data, data_imp, data_pca_multilevel,
    pca_withna, pca_impute, mdist, argv, file=argv$rdata
  )

  # partial least squares discriminant analysis
  if (argv$plsdacomp > 0) {
    if (!is.na(pch)) {
      data_plsda = plsda_classify(
        input_data, classes, pch, title=names, argv$plsdacomp
      )
    } else {
      data_plsda = plsda_classify(
        input_data, classes, pch=NA, title=names, argv$plsdacomp
      )
    }
  } else { data_plsda = NA }

  save(classes, data, input_data, data_pca_multilevel, data_plsda, pca_withna,
    pca_impute, mdist, argv, file=argv$rdata
  )

  # sparse partial least squares discriminant analysis
  if (argv$splsdacomp > 0) {
    if (!is.na(pch)) {
      print("Tuning splsda components and selected variables")
      if (is.na(argv$splsda_keepx)) {
        splsda_keepx = c(1,2,3)
        splsda_ncomp = length(splsda_keepx)
      } else {
        splsda_keepx = lapply(strsplit(argv$splsda_keepx, ","), as.integer)[[1]]
        splsda_ncomp = length(splsda_keepx)
      }

      print("sPLSDA keepX:")
      print(splsda_keepx)
      print("sPLSDA ncomp:")
      print(splsda_ncomp)

      tuned = splsda_tune(input_data, classes, names, data.frame(pch),
        ncomp=splsda_ncomp, nrepeat=10, logratio="none",
        test_keepX=splsda_keepx, validation="loo", folds=10, dist=argv$mdist,
        cpus=argv$ncpus, progressBar=TRUE)

      splsda_keepx = lapply(tuned, `[`, "choice.keepX")
      splsda_ncomp = lapply(tuned, `[`, "choice.ncomp")

      print("Tuned splsda to use number of components:")
      splsda_ncomp = lapply(splsda_ncomp, `[`, "ncomp")
      splsda_ncomp = unlist(splsda_ncomp, recursive = FALSE)
      names(splsda_ncomp) = names
      print(splsda_ncomp)

      print("Tuned the number of variables selected on each component to:")
      print(splsda_keepx)
      splsda_keepx = unlist(splsda_keepx, recursive = FALSE)
      names(splsda_keepx) = names
      print(splsda_keepx)

      data_splsda = splsda_classify(
        data_imp, classes, pch, title=names, splsda_ncomp, splsda_keepx
      )
    }
  } else {
    data_splsda = NA
    tuned = NA
  }

  save(classes, data, data_imp, data_pca_multilevel, data_plsda, data_splsda,
    tuned, pca_withna, pca_impute, mdist, argv, file=argv$rdata
  )
  q()
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
