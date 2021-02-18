#!/usr/bin/Rscript
# regenerate plots for this manuscript
library(mixOmics)
library(reshape2)

make_indiv <- function(data, classes, title, dist="centroids.dist", 
		       pch=NA, bg=FALSE) {
  # take a mixomics data object and make plots
  if (bg == TRUE) {
    bg <- background.predict(data, comp.predicted=2, dist=dist)
  } else { bg <- NULL }
  if (is.na(pch)) {
    pch <- seq(unique(classes))
  }
  data <- plotIndiv(
    data, ind.names=FALSE, group=classes, legend=TRUE,
    pch=pch, title=paste(title, "1/2"), comp=c(1,2), ellipse=TRUE,
    background=bg
  )
  return(data)
}

make_cim <- function(data, classes=NA, title="", pch=NA, 
		     show_cols=FALSE) {
  # take a mixomics data object and make clustered image maps
  colours_class = color.mixo(
      1:length(unique(classes))
    )[as.numeric(as.factor(classes))]
  if (!is.na(pch)) {
    colours_pch = color.mixo(
      1:length(unique(pch))
    )[as.numeric(as.factor(pch))]
    colours_cim = cbind(colours_class, colours_pch)
  } else {colours_cim = data.frame(colours_class)}

  is_diablo <- c("block.splsda", "block.spls", "sgccda", "sgcca", "DA")
  if (all(class(data) == is_diablo)) {
    data <- cimDiablo(
      data, title=title, col.names=show_cols, margin=c(6,12),
      size.legend=0.5
    )
  } else {
    data <- cim(data, title=title, row.sideColors=colours_cim,
      legend=list(title="Status"), col.names=show_cols,
      margin=c(6,12)
    )
  }
  return(data)
}

make_loadings <- function(data, comp=1, title="", short=6, 
			  force_short=FALSE) {
  # take a mixomics data object and plot top N loadings
  name_var <- data$names$colnames$X
  if (!short == 0) {
    name_var <-  make.names(sapply(
      name_var, strtrim, short, USE.NAMES=FALSE
      ), unique=TRUE)
  }

  is_diablo <- c("block.splsda", "block.spls", "sgccda", "sgcca", "DA")
  if (all(class(data) == is_diablo)) {
    name_var <- colnames(data)
    size_name <- 0.001
  } else { 
    size_name <- 1 
  }

  if (force_short == TRUE) {
    data_vis_names = replace_names_(data, trim=6)
    data_vis = data_vis_names$data_vis
    data <- data_vis
  }

  data <- plotLoadings(
      data, contrib="max", comp=comp, max.name.length=6,
      method='median', ndisplay=20, name.var=name_var, 
      size_name=size_name, size.legend=0.6, title=title,
      max.name.length=6, 
      )
  return(data)
}

make_error <- function(data, title="", cpus=6, stability=FALSE) {
  # take a mixomics data object and calculate error rates
  # and variable stability (diablo is hard-coded to show)
  is_diablo <- c("block.splsda", "block.spls", "sgccda", "sgcca", "DA")
  is_diablo <- all(class(data) == is_diablo)
  data <- perf(
    object=data, validation="loo", nrepeats=10, auc=TRUE, cpus=cpus,
    progressBar=TRUE
  )
  data_error <- plot(data, main=paste("Error rate", title), 
       col=color.mixo(5:7), sd=TRUE)
  if (stability == TRUE && is_diablo == FALSE) {
    data_stability <- plot(data$features$stable[[1]], 
      type="h", main=paste("Stability Comp 1", title), las=2,
      ylab="Stability", xlab="Features", xaxt='n'
    )
  } else {data_stability <- NA}
  if (stability == TRUE && is_diablo == TRUE) {
    blocks <- names(data$features$stable$nrep1)
    for (i in blocks) {
      j <- data$features$stable$nrep1[[i]]$comp1
      plot(j, type="h", las=2, ylab="Stability", xlab="Features",
        main=paste(i, "Stability Comp 1", title), xaxt='n'
      )
    }
  } else {data_stability <- NA}
  return(list(data_error, data_stability))
}

make_tuned <- function(data, title="") {
  # take tuned mixomics object and plot error rates
  data <- plot(data, main=paste("Error rate", title))
  return(data)
}

make_auroc <- function(data, title="", pch=NA) {
  # take tuned mixomics object and plot roc curve
  if (is.na(pch)) {pch <- NULL}
  data <- auroc(data, title=paste("ROC comp 1", title),
  		multilevel=pch)
  return(data)
}

make_all <- function(data, title, classes, outfile_path, 
		     pch=NULL, cpus=6, dist="centroids.dist", 
		     show_cols=FALSE, tuned=NA, bg=FALSE,
		     stability=FALSE, force_short=FALSE, auc=TRUE) {
  # aggregate all functions into master plotting function
  print(paste("Making plots for:", title))
  pdf(outfile_path)
  make_indiv(data=data, classes=classes, dist=dist,
             title=title, pch=pch, bg=bg)
  make_cim(data=data, classes=classes, show_cols=show_cols,
           title=title, pch=pch)
  make_loadings(data=data, title=title, force_short=force_short)
  make_error(data=data, title=title, stability=stability)
  if (!is.na(tuned)) {make_tuned(data=tuned, title=title)}
  if (auc == TRUE) {make_auroc(data=data, title=title, pch=pch)}
  dev.off()
  print(paste("Plots written to:", outfile_path))
}

# replace names in all associated columns for visualisation only
replace_names_ <- function(data, trim=16) {
  all_names <- head(data$names$colnames, n=-1)
  split_ <- function(block, all_names) {
    all_names <- gsub("__FEATUREID", "", all_names[[block]])
    all_names <- gsub(paste("_", block, sep=""), "", all_names)
    return(all_names)
  }
  splitted <- mapply(
    function(x) split_(x, all_names), names(all_names), SIMPLIFY=FALSE
  )
  truncate_ <- function(names, trim) {
    ifelse(nchar(names) > trim, paste0(strtrim(names, trim), ''), names)
  }
  truncated <- lapply(splitted, truncate_, trim)
  # make a copy, dont want to overwrite
  data_vis <- data
  # map truncated values to all locations
  for (i in names(all_names)) {
    row.names(data_vis$loadings[[i]]) <- make.unique(
      sapply(truncated[[i]], toString), sep="__"
    )
    data_vis$names$colnames[[i]] <- make.unique(
      sapply(truncated[[i]], toString), sep="__"
    )
    colnames(data_vis$X[[i]]) <- make.unique(
      sapply(truncated[[i]], toString), sep="__"
    )
  }
  return(list(data_vis=data_vis, truncated=truncated))
}

show_na_prop <- function(data_na, name) {
  # first, sum the number of missing values per variable
  sum_na_per_var <- apply(data_na, 2, function(x) {sum(is.na(x))})
  # show proportion of NA values across all samples (y) for a variable (x)
  plot(sum_na_per_var/nrow(data_na), type='h', xlab='variable index',
    ylab='NA proportion (across all samples for a variable)',
    main=paste(name, 'NA rate per variable on unfiltered data'))
}

cor_imputed_unimputed = function(pca_withna, pca_impute, names) {
  # plots a heatmap of correlations: -> df, df, name
  x <- pca_withna
  y <- pca_impute
  z <- names
  print("Plotting correlation between unimputed and imputed components")
  print(ggplot(melt(cor(x$variates$X, y$variates$X)),
    aes(Var1, Var2, fill=value)) +
    ggtitle(paste(z, "Correlation between imputed and unimputed data")) +
    geom_tile() +
    scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0,
      limit=c(-1,1), space="Lab", name="Pearson\nCorrelation") +
    theme_minimal()
  )
}

make_case_study_1 <- function() {
  infile_path_1 <- "../results/RData.RData"
  outfile_dir <- "../results/manuscript_figures/"
  
  # case study 1
  load(infile_path_1)

  # convert to mixOmics data objects for further operations
  # (for compatibility with the original iteration of code)
  proteome_plsda <- data_plsda[,1]
  class(proteome_plsda) <- c("mixo_plsda", "mixo_pls", "DA")

  translatome_plsda <- data_plsda[,2]
  class(translatome_plsda) <- c("mixo_plsda", "mixo_pls", "DA")

  proteome_splsda <- data_splsda[,1]
  class(proteome_splsda) <- c("mixo_splsda", "mixo_spls", "DA")

  translatome_splsda <- data_splsda[,2]
  class(translatome_splsda) <- c("mixo_splsda", "mixo_spls", "DA")

  make_all(
    outfile_path=paste(outfile_dir, "case_1_plsda_prot.pdf", sep=""),
    data=proteome_plsda, title="Proteome PLSDA", classes=classes,
    pch=pch, cpus=6, dist="centroids.dist", show_cols=FALSE, 
    tuned=NA, bg=TRUE, stability=FALSE
  )

  make_all(
    outfile_path=paste(outfile_dir, "case_1_plsda_tran.pdf", sep=""),
    data=translatome_plsda, title="Translatome PLSDA", classes=classes,
    pch=pch, cpus=6, dist="centroids.dist", show_cols=FALSE, 
    tuned=NA, bg=TRUE, stability=FALSE
  )

  make_all(
    outfile_path=paste(outfile_dir, "case_1_splsda_prot.pdf", sep=""),
    data=proteome_splsda, title="Proteome sPLSDA", classes=classes,
    pch=pch, cpus=6, dist="centroids.dist", show_cols=FALSE, 
    tuned=tuned_splsda$proteome, bg=TRUE, stability=TRUE
  )

  make_all(
    outfile_path=paste(outfile_dir, "case_1_splsda_tran.pdf", sep=""),
    data=translatome_splsda, title="Translatome sPLSDA", classes=classes,
    pch=pch, cpus=6, dist="centroids.dist", show_cols=FALSE, 
    tuned=tuned_splsda$translatome, bg=TRUE, stability=TRUE
  )

  make_all(
    outfile_path=paste(outfile_dir, "case_1_diablo.pdf", sep=""),
    data=diablo, title="Multiblock sPLSDA", classes=classes,
    pch=pch, cpus=6, dist="mahalanobis.dist", show_cols=FALSE,
    tuned=NA, bg=FALSE, force_short=TRUE, stability=TRUE
  )

  # clear environment to avoid clashes with case study 2
  rm(list=ls())
}

make_case_study_1_extra <- function() {
  # case study 1 extra plots
  infile_path_1 <- "../results/reformatted_output/RData.RData"
  unimputed_prot_path <- "../data/diablo_proteome.tsv"
  unimputed_prot <- read.table(unimputed_prot_path, sep="\t", header=TRUE, row.names=1)
  unimputed_tran_path <- "../data/diablo_translatome.tsv"
  unimputed_tran <- read.table(unimputed_tran_path, sep="\t", header=TRUE, row.names=1)
  outfile_dir <- "../results/manuscript_figures/"

  # case study 1
  load(infile_path_1)
  
  pdf(paste(outfile_dir, "case_1_extra_plots.pdf"))
  
  unimputed_prot[unimputed_prot == 0] <- NA
  unimputed_tran[unimputed_tran == 0] <- NA
  na_prop_prot <- show_na_prop(
    unimputed_prot, "Missing values in original proteome data"
  )
  na_prop_tran <- show_na_prop(
    unimputed_tran, "Missing values in original translatome data"
  )
  pca_unimputed_prot <- pca(unimputed_prot, ncomp=10)
  pca_unimputed_tran <- pca(unimputed_tran, ncomp=10)
  pca_imputed_prot <- pca(data$proteome, ncomp=10)
  pca_imputed_tran <- pca(data$translatome, ncomp=10)
  
  prot <- "Proteome correlation imputed vs unimputed data"
  tran <- "Translatome correlation imputed vs unimputed data"
  cor_imputed_unimputed(
    pca_imputed_prot, pca_unimputed_prot, prot
  )
  cor_imputed_unimputed(
    pca_imputed_tran, pca_unimputed_tran, tran
  )
  dev.off()

  # clear environment
  rm(list=ls())
}

make_case_study_2 <- function() {
  
  infile_path_2 <- "../results/MSV000085703/reformatted_output/RData.RData"
  outfile_dir <- "../results/manuscript_figures/"
  
  # case study 2
  load(infile_path_2)
  
  lipi_plsda_data <- data_plsda$lipidome$data_plsda
  meta_plsda_data <- data_plsda$metabolome$data_plsda
  prot_plsda_data <- data_plsda$proteome$data_plsda
  tran_plsda_data <- data_plsda$transcriptome$data_plsda
  omic_plsda_data <- list(
    lipi_plsda_data, meta_plsda_data, prot_plsda_data, tran_plsda_data
  )

  lipi_plsda_perf <- data_plsda$lipidome$perf_plsda
  meta_plsda_perf <- data_plsda$metabolome$perf_plsda
  prot_plsda_perf <- data_plsda$proteome$perf_plsda
  tran_plsda_perf <- data_plsda$transcriptome$perf_plsda
  omic_plsda_perf <- list(
    lipi_plsda_perf, meta_plsda_perf, prot_plsda_perf, tran_plsda_perf
  )

  lipi_splsda_data <- data_splsda$lipidome$data_splsda
  meta_splsda_data <- data_splsda$metabolome$data_splsda
  prot_splsda_data <- data_splsda$proteome$data_splsda
  tran_splsda_data <- data_splsda$transcriptome$data_splsda
  omic_splsda_data <- list(
    lipi_splsda_data, meta_splsda_data, prot_splsda_data, tran_splsda_data
  )

  lipi_splsda_perf <- data_splsda$lipidome$perf_splsda
  meta_splsda_perf <- data_splsda$metabolome$perf_splsda
  prot_splsda_perf <- data_splsda$proteome$perf_splsda
  tran_splsda_perf <- data_splsda$transcriptome$perf_splsda
  omic_splsda_perf <- list(
    lipi_splsda_perf, meta_splsda_perf, prot_splsda_perf, tran_splsda_perf
  )

  blocks_plsda <- names(data_plsda)
  mapply(function(x, y) make_all(
    outfile_path=paste(outfile_dir, "case_2_plsda_", y, ".pdf", sep=""),
    data=x, title=paste(y, "PLSDA"), classes=classes, pch=pch,
    cpus=6, dist="centroids.dist", show_cols=FALSE, bg=TRUE,
    tuned=NA, stability=FALSE
    ),
    omic_plsda_data, blocks_plsda, SIMPLIFY=FALSE
  )

  blocks_splsda <- names(data_splsda)
  mapply(function(x, y) make_all(
    outfile_path=paste(outfile_dir, "case_2_splsda_", y, ".pdf", sep=""),
    data=x, title=paste(y, "sPLSDA"), classes=classes, pch=pch,
    cpus=6, dist="centroids.dist", show_cols=FALSE, bg=TRUE,
    tuned=NA, stability=TRUE
    ),
    omic_splsda_data, blocks_splsda, SIMPLIFY=FALSE
  )

  make_all(
    outfile_path=paste(outfile_dir, "case_2_diablo.pdf", sep=""),
    data=diablo, title="Multiblock sPLSDA", classes=classes,
    pch=pch, cpus=6, dist="mahalanobis.dist", show_cols=FALSE,
    tuned=NA, bg=FALSE, force_short=TRUE, stability=TRUE, auc=FALSE
  )

  # clear environment
  rm(list=ls())
}

make_case_study_2_extra <- function() {
  # case study 2 extra plots
  infile_path_2 <- "../results/MSV000085703/reformatted_output/RData.RData"
  unimputed_path <- "../data/MSV000085703/hfd45/data_transcriptomics.tsv"
  unimputed <- read.table(unimputed_path, sep="\t", header=TRUE, row.names=1)
  outfile_dir <- "../results/manuscript_figures/"

  # case study 2
  load(infile_path_2)

  pdf(paste(outfile_dir, "case_2_extra_plots.pdf", sep=""))

  unimputed[unimputed == 0] <- NA
  na_prop <- show_na_prop(
    unimputed, "Missing values in original transcriptome data"
  )
  pca_imputed <- pca(data$transcriptome, ncomp=10)
  pca_unimputed <- pca(unimputed, ncomp=10)
  
  name <- "Transcriptome correlation imputed vs unimputed data"
  cor_imputed_unimputed(pca_imputed, pca_unimputed, name)
  dev.off()

  # clear environment
  rm(list=ls())
}

main <- function() {
  # case study 1 was generated with an old version of the code
  # some data structure change but underlying info is the same
  make_case_study_1()
  make_case_study_1_extra()  

  # case study 2 was generated with an up to date code version
  make_case_study_2()
  make_case_study_2_extra()
}

main()
