#!/usr/bin/Rscript
# combine multi-omics data
library(mixOmics)
library(argparser, quietly=TRUE)
library(ggplot2)
library(multiomics)
library(parallel)
library(reshape2)

parse_argv <- function() {
  library(argparser, quietly=TRUE)
  p <- argparser::arg_parser(
    "Assess the quality of individual omics data modality blocks, \
    Run DIABLO on multi-omics data. Take tsv file of classes, tsv files of omics \
    data (at least 2) and identify correlation between features. \
    For more information please refer to: \
      our pipeline documentation at https://gitlab.com/tyagilab/sars-cov-2/ \
      mixOmics and case studies at http://mixomics.org/mixdiablo/case-study-tcga/"
  )
  # Add command line arguments
  p <- argparser::add_argument(
    p, "--json", type="character", default=NA,
    help="pass args as json file instead of command line args (overrides args!)"
  )
  p <- argparser::add_argument(
    p, "--classes", help="sample information", type="character"
  )
  p <- argparser::add_argument(
    p, "--classes_secondary", type="character", default=NA,
    help="secondary sample information eg individual (same format as classes)"
  )
  p <- argparser::add_argument(
    p, "--dropna_classes", type="bool", default=FALSE,
    help="where all replicates for >= 1 class are NA, drop that feature. \
    If both --dropna_classes and --dropna_prop are enabled, perform \
    --dropna_classes and then --dropna_prop."
  )
  p <- argparser::add_argument(
    p, "--dropna_prop", type="integer", default=0.6,
    help="drop feature that does not meet NA proportion threshold, eg if 0.3, \
    drop a feature if >= 0.3 of values are NA. If both --dropna_classes and \
    --dropna_prop are enabled, perform --dropna_classes and then --dropna_prop."
  )
  p <- argparser::add_argument(
    p, "--data", type="character", nargs=Inf, default=NA,
    help="paths to omics data (must match order of names in --data_names)."
  )
  p <- argparser::add_argument(
    p, "--data_names", type="character", nargs=Inf, default=NA,
    help="names of the individual omics data blocks (must match order of data in --data)."
  )
  p <- argparser::add_argument(
    p, "--force_unique", type="bool", default=TRUE,
    help="force values to be unique in each omics data block (default TRUE)"
  )
  p <- argparser::add_argument(
    p, "--mappings", type="character", nargs=Inf,
    help="path to map file of feature id to name (must be same order as data!)"
  )
  p <- argparser::add_argument(
    p, "--ncpus", help="number of cpus", type="integer", default=parallel::detectCores()-2
  )
  p <- argparser::add_argument(
    p, "--diablocomp", type="integer", default=0,
    help="number of diablo components (set manually if you get inference error). Default to number of classes - 1."
  )
  p <- argparser::add_argument(
    p, "--linkage", type="integer", default=0.1,
    help="degree of blocks linkage 0 (discrimination) < x < 1 (correlation)"
  )
  p <- argparser::add_argument(
    p, "--diablo_keepx", type="vector", default=NA, nargs="+",
    help="variables to keep for diablo"
  )
  p <- argparser::add_argument(
    p, "--icomp", type="integer", default=0,
    help="component number for imputing (set 0 for no imputation)"
  )
  p <- argparser::add_argument(
    p, "--zero_as_na", default=TRUE,
    help="treat zero values as missing values for imputation (DEFAULT: TRUE)"
  )
  p <- argparser::add_argument(
    p, "--replace_missing", default=TRUE,
    help="replace missing values only in imputation, else replace all values (DEFAULT: TRUE)"
  )
  p <- argparser::add_argument(
    p, "--pcomp", type="integer", default=5,
    help="number of principal components (defaults to 5)"
  )
  p <- argparser::add_argument(
    p, "--plsdacomp", type="integer", default=0,
    help="number of components for plsda (defaults to number of samples)"
  )
  p <- argparser::add_argument(
    p, "--splsdacomp", type="integer", default=0,
    help="number of components for splsda (defaults to number of samples)"
  )
  p <- argparser::add_argument(
    p, "--splsda_keepx", type="vector", default=NA, nargs="+",
    help="variables to keep for splsda"
  )
  p <- argparser::add_argument(
    p, "--dist_plsda", type="character", default="centroids.dist",
    help="plsda distance metric. [max.dist, centroids.dist, mahalanobis.dist]"
  )
  p <- argparser::add_argument(
    p, "--dist_splsda", type="character", default="centroids.dist",
    help="splsda distance metric [max.dist, centroids.dist, mahalanobis.dist]"
  )
  p <- argparser::add_argument(
    p, "--dist_diablo", type="character", default="centroids.dist",
    help="diablo distance metric [max.dist, centroids.dist, mahalanobis.dist]"
  )
  p <- argparser::add_argument(
    p, "--cross_val", type="character", default="loo",
    help="cross-validation method ['Mfold', 'loo'], defaults to loo"
  )
  p <- argparser::add_argument(
    p, "--cross_val_nrepeat", type="integer", default=10,
    help="cross-validation, repeat this many times"
  )
  p <- argparser::add_argument(
    p, "--cross_val_folds", type="integer", default=10,
    help="folds in cross-validation, applicable to Mfold only"
  )
  p <- argparser::add_argument(
    p, "--contrib", type="character", default="max",
    help="contribution type for plotting loadings of s/PLSDA/DIABLO [max|min]"
  )
  p <- argparser::add_argument(
    p, "--corr_cutoff", type="double", default=0.95,
    help="correlation cutoff for displaying lines on circos plot"
  )
  p <- argparser::add_argument(
    p, "--low_var", flag=TRUE,
    help="enable low variance mode (if you have low variance or sparse data)"
  )
  p <- argparser::add_argument(
    p, "--tune_off", flag=TRUE,
    help="turn off tuning (if you already know the values or for debugging)"
  )
  p <- argparser::add_argument(
    p, "--outfile_dir", type="character", default="./",
    help="write args, R plots and RData here (will overwrite existing!)"
  )
  p <- argparser::add_argument(
    p, "--rdata", type="character", default="./data.RData",
    help="write RData object here, has (classes, data, diablo, mdist)"
  )
  p <- argparser::add_argument(
    p, "--plot", type="character", default="./Rplots.pdf",
    help="write R plots here (will overwrite existing!)"
  )
  p <- argparser::add_argument(
    p, "--args", type="character", default="Rscript.sh",
    help="command line options for script are saved here as a shell file"
  )
  p <- argparser::add_argument(
    p, "--mini_run", flag=TRUE,
    help="whether to run on first 100 features (max) from each dataset"
  )
  # Parse the command line arguments
  argv <- argparser::parse_args(p)

  # Do work based on the passed arguments
  return(argv)
}

write_args <- function(args, argpath) {
  args <- as.data.frame(stack(args))
  args <- args[4:dim(args)[1],]
  args <- args[,c(2,1)]
  args <- aggregate(args$values, list(args$ind), paste, collapse=" ")
  colnames(args) <- c("ind", "values")
  args <- paste("--", paste(args$ind, args$values), " \\", sep="")
  script_name <- paste("Rscript", sub(".*=", "", commandArgs()[4]), "\\")
  for (i in grep(",", args)) { args[i] = gsub(",", " ", args[i]) }
  last <- paste(unlist(strsplit(args[length(args)], " "))[2], " \\")
  args <- head(args, -1)
  args[length(args)] <- substr(args[length(args)],1,nchar(args[length(args)])-2)
  args <- c(script_name, paste("  ", last), paste("  ", args))
  write(args, sep="\n", file=argpath)
}

main <- function() {
  argv <- parse_argv()

  if (!is.na(argv$json)) {
    print("Json file passed to pipeline, override all other command input!")
    argv <- rjson::fromJSON(file=argv$json)
  }

  print("Creating output files directory (will overwrite existing data!)")
  outdir <- argv$outfile_dir
  dir.create(file.path(outdir), recursive = TRUE, showWarnings = FALSE)

  print("Writing command line arguments to:")
  argpath <- paste(outdir, argv$args, sep="/")
  print(argpath)
  write_args(argv, argpath)

  rdata <- paste(outdir, argv$rdata, sep="/")
  plot <- paste(outdir, argv$plot, sep="/")
  contrib <- argv$contrib

  # print some diagnostics for debugging
  print("Available cpus:")
  print(parallel::detectCores())
  print("Using cpus (change with --ncpus):")
  print(argv$ncpus)
  print("Degree of linkage for omics blocks:")
  linkage <- argv$linkage
  print(linkage)
  print("Use low variance mode (for sparse or low variance data):")
  low_var <- argv$low_var
  dist_plsda <- argv$dist_plsda
  dist_splsda <- argv$dist_splsda
  dist_diablo <- argv$dist_diablo
  corr_cutoff <- argv$corr_cutoff
  print("Perform tuning of sPLSDA and DIABLO:")
  tune_off <- argv$tune_off
  print(!tune_off)
  print("Cross-validation method:")
  print(argv$cross_val)
  print("Cross-validation nrepeat:")
  print(argv$cross_val_nrepeat)
  print("Cross-validation folds:")
  print(argv$cross_val_folds)
  print("Distance measure (PLSDA):")
  print(dist_plsda)
  print("Distance measure (sPLSDA):")
  print(dist_splsda)
  print("Distance measure (DIABLO):")
  print(dist_diablo)
  print("Display correlation cutoff:")
  print(corr_cutoff)
  options(warn=1)

  paths <- argv$data
  print("Paths to data:")
  print(paths)
  print("Parsing classes")
  classes <- parse_classes(argv$classes)

  if (!is.na(argv$classes_secondary)) {
    print("Parsing secondary classes")
    pch <- parse_classes(argv$classes_secondary)
  } else {
    pch <- NA
    tuned_splsda <- NA
    data_splsda <- NA
  }

  data_names <- argv$data_names
  print("Omics data types")
  print(data_names)

  if (length(data_names) != length(paths) | any(duplicated(data_names)) | any(duplicated(paths))) {
    stop("Data names must correspond to paths and be unique!")
    print("Data names:")
    print(data_names)
    print("Data paths:")
    print(paths)
    quit(save="no", status=1)
  }

  # initialise plots
  print(paste("Saving plots to (overwriting existing):", plot))
  pdf(plot)

  # load data and drop features / columns with all NA
  data <- lapply(paths, parse_data, missing_as=NA, rmna=TRUE)
  names(data) <- data_names

  ## mini run
  if (isTRUE(argv$mini_run))
  {
    print('Performing a mini run as --mini_run flag is used...')
    data <- lapply(data, function(x) {
      minirun_ncol <- min(300, ncol(x))
      x[,seq_len(minirun_ncol)]
    })
  }


  # show proportion of NA values in unfiltered data
  for (i in length(data))
  {
    x <- data[[i]]
    y <- data_names[[i]]
    show_na_prop(x, y)
  }
  # drop features / columns where >= 1 class is not represented
  if (argv$dropna_classes == TRUE) {
    data <- lapply(data, remove_na_class, classes)
    save(list = ls(all.names = TRUE), file=rdata)
  }

  # drop features / columns >= a threshold of NA values
  if (argv$dropna_prop > 0) {
    data <- remove_na_prop(
      data, na_prop=argv$dropna_prop, argv$zero_as_na
    )
    save(list = ls(all.names = TRUE), file=rdata)
  }

  if (length(argv$mappings) > 1) {
    print("Using mappings from files (order must be identical to data!):")
    print(argv$mappings)
    mappings <- argv$mappings
    mapped <- lapply(mappings, parse_mappings)
    names(mapped) <- data_names
    data <- mapply(function(x, y) remap_data(x, y), data, mapped)
    names(data) <- data_names
  } else {
    print("Not remapping new feature names to existing, will use original.")
    mappings <- NA
  }

  if (argv$force_unique == TRUE) {
    diablo_input <- force_unique_blocks(data)
  }

  # check dimensions
  print("Data dimensions:")
  dimensions <- lapply(data, dim)
  print(dimensions)

  design <- create_design(data, linkage)

  # check classes
  print(summary(classes))
  print("Y (classes):")
  print(table(classes))
  print("Design:")
  print(design)

  # count missing data after all filtering
  missing <- lapply(data, count_missing)
  pca_withna <- plot_pca_single(
    data, classes, pch=pch, ncomp=argv$pcomp,
    title=paste("No Impute n_PCs =", argv$pcomp, "\n")
  )
  save(list = ls(all.names = TRUE), file=rdata)

  # impute data if components given
  # refer to http://mixomics.org/methods/missing-values/
  print(argv$icomp)
  if (argv$icomp > 0) {
    print("Impute components set, imputing NA values (set -i 0 to disable)")
    data_imp <- impute_missing(data, rep(argv$icomp, length(data)), outdir)

    if (argv$replace_missing == TRUE) {
      data_imp <- replace_missing(data, data_imp)
    }

    pca_impute <- plot_pca_single(
      data_imp, classes, pch=pch, ncomp=argv$pcomp,
      title=paste("Imputed n_PCs =", argv$pcomp, "IC:", argv$icomp, "\n")
    )
    heatmaps <- cor_imputed_unimputed(pca_withna, pca_impute, data_names)
  } else {
    print("Impute components unset, not imputing NA (set -i > 0 to enable)")
    data_imp <- NA
    pca_impute <- NA
  }
  save(list = ls(all.names = TRUE), file=rdata)

  # multilevel decomposition if secondary variables are specified
  # refer to http://mixomics.org/case-studies/multilevel-vac18/
  # multilevel pca
  if (!is.na(data_imp)) {
    input_data <- data_imp
  } else {
    input_data <- data
  }

  if (!is.na(argv$classes_secondary)) {
    data_pca_multilevel <- plot_pca_multilevel(
      input_data, classes, pch=pch, ncomp=argv$pcomp,
      title=paste("No Impute. n_PCs =", argv$pcomp, "\n")
    )
    if (exists("data_imp")) {
      data_pca_multilevel <- plot_pca_multilevel(
        input_data, classes, pch=pch, ncomp=argv$pcomp,
        title=paste("Imputed. n_PCs =", argv$pcomp, "IC:", argv$icomp, "\n")
      )
    }
  } else { data_pca_multilevel <- NA }
  save(list = ls(all.names = TRUE), file=rdata)

  # partial least squares discriminant analysis
  if (argv$plsdacomp > 0) {
    if (length(pch) > 1) {
      data_plsda <- classify_plsda(input_data, classes, pch, title=data_names,
        argv$plsdacomp, contrib, outdir, mappings, dist_splsda, bg=TRUE,
        validation=argv$cross_val, folds=argv$cross_val_folds,
        nrepeat=argv$cross_val_nrepeat, near_zero_var=low_var
      )
    } else {
      data_plsda <- classify_plsda(input_data, classes, pch=NA, title=data_names,
        argv$plsdacomp, contrib, outdir, mappings, dist_splsda, bg=TRUE,
        validation=argv$cross_val, folds=argv$cross_val_folds,
        nrepeat=argv$cross_val_nrepeat, near_zero_var=low_var
      )
    }
  } else { data_plsda <- NA }

  save(list = ls(all.names = TRUE), file=rdata)

  # sparse partial least squares discriminant analysis
  if (argv$splsdacomp > 0) {
      if (!is.na(argv$splsda_keepx)) {
        splsda_keepx <- lapply(strsplit(argv$splsda_keepx, ","), as.integer)[[1]]
        # splsda_ncomp = length(splsda_keepx)
      }

      splsda_ncomp <- argv$splsdacomp
      print("sPLSDA keepX:")
      print(splsda_keepx)
      print("sPLSDA ncomp:")
      print(splsda_ncomp)

      if (!tune_off) {
        print("Tuning splsda components and selected variables")
        if (length(pch) > 1) {
          tuned_splsda <- tune_splsda(input_data, classes, data_names,
            data.frame(pch),
            ncomp=splsda_ncomp, nrepeat=argv$cross_val_nrepeat, logratio="none",
            test_keepX=splsda_keepx, validation=argv$cross_val,
            folds=argv$cross_val_folds, dist=dist_splsda, cpus=argv$ncpus,
            progressBar=TRUE, near_zero_var=low_var)
        } else {
          tuned_splsda <- tune_splsda(input_data, classes, data_names,
            NULL,
            ncomp=splsda_ncomp, nrepeat=argv$cross_val_nrepeat, logratio="none",
            test_keepX=splsda_keepx, validation=argv$cross_val,
            folds=argv$cross_val_folds, dist=dist_splsda, cpus=argv$ncpus,
            progressBar=TRUE, near_zero_var=low_var)
        }
        splsda_keepx <- lapply(tuned_splsda, `[`, "choice.keepX")
        splsda_ncomp <- lapply(tuned_splsda, `[`, "choice.ncomp")

        print("Tuned splsda to use number of components:")
        splsda_ncomp <- lapply(splsda_ncomp, `[`, "ncomp")
        splsda_ncomp <- unlist(splsda_ncomp, recursive = FALSE)
        names(splsda_ncomp) <- data_names
        print(splsda_ncomp)

        print("Tuned the number of variables selected on each component to:")
        print(splsda_keepx)
        splsda_keepx <- unlist(splsda_keepx, recursive = FALSE)
        names(splsda_keepx) <- data_names
      } else {
        print("No sPLSDA tuning performed!")
      }

      if (length(pch) > 1) {
        data_splsda <- classify_splsda(
          input_data, classes, pch, title=data_names, splsda_ncomp,
          splsda_keepx, contrib, outdir, mappings, data_splsda, bg=TRUE,
          near_zero_var=low_var
        )
      } else {
        data_splsda <- classify_splsda(
          input_data, classes, pch=NA, title=data_names, splsda_ncomp,
          splsda_keepx, contrib, outdir, mappings, data_splsda, bg=TRUE,
          near_zero_var=low_var
        )
      }
  } else {
    data_splsda <- NA
    tuned_splsda <- NA
  }

  save(list = ls(all.names = TRUE), file=rdata)

  # NOTE: if you get tuning errors, set dcomp manually with --dcomp N
  if (!tune_off) {
    tuned_diablo <- tune_diablo_ncomp(
      data, classes, design, argv$diablocomp, cpus=argv$ncpus,
      near_zero_var=low_var
    )
    perf_diablo <- tuned_diablo
    print("Parameters with lowest error rate:")
    tuned_diablo <- tuned_diablo$choice.ncomp$WeightedVote["Overall.BER",]
    diablo_ncomp <- tuned_diablo[which.max(tuned_diablo)]
  } else {
    print("No DIABLO tuning (number of components) performed!")
    diablo_ncomp <- argv$diablocomp
    tuned_diablo <- NA
    perf_diablo <- NA
  }
  print("Number of components:")
  print(diablo_ncomp)

  # remove invariant columns
  # data = lapply(data, remove_novar)
  save(list = ls(all.names = TRUE), file=rdata)

  # block-wise splsda doesnt do internal multilevel decomposition
  if (length(pch) > 1) {
    diablo_input <- lapply(input_data, withinVariation, design=data.frame(pch))
  } else {
    diablo_input <- input_data
  }

  print("Making feature names unique across all blocks...")
  if (argv$force_unique == TRUE) {
    diablo_input <- force_unique_blocks(data)
  }
  save(list = ls(all.names = TRUE), file=rdata)

  # tune diablo parameters and run diablo
  diablo_keepx <- lapply(strsplit(argv$diablo_keepx, ","), as.integer)[[1]]

  if (!tune_off) {
    diablo_keepx <- tune_diablo_keepx(diablo_input, classes, diablo_ncomp,
      design, diablo_keepx, cpus=argv$ncpus, dist=dist_diablo, progressBar=TRUE,
      validation=argv$cross_val, folds=argv$cross_val_folds,
      nrepeat=argv$cross_val_nrepeat, near_zero_var=low_var
    )
    print("Diablo keepx:")
    print(diablo_keepx)
    diablo <- run_diablo(
      diablo_input, classes, diablo_ncomp, design, diablo_keepx, low_var
    )
  } else {
    print("No DIABLO tuning (keepX) performed!")
    diablo_keepx <- rep(diablo_keepx, length(data_names))
    diablo_keepx <- split(
      diablo_keepx, sort(rep_len(1:length(data_names), length(diablo_keepx)))
    )
    names(diablo_keepx) <- data_names
    print("Diablo keepx:")
    print(diablo_keepx)
    diablo <- run_diablo(
      diablo_input, classes, diablo_ncomp, design, diablo_keepx, low_var
    )
  }

  print("Diablo design:")
  print(diablo$design)
  plot_diablo(diablo, diablo_ncomp, outdir, data_names, "keepx", corr_cutoff)
  # assess_performance(diablo, dist=dist_diablo, diablo_ncomp)

  # save RData object for future reference
  save(list = ls(all.names = TRUE), file=rdata)
  dev.off()
}

main()
