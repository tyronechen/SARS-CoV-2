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
    ylab='missing values proportion (across all samples for a variable)',
    main=paste(name, 'missing values per variable on unfiltered data'))
}

remove_na_class <- function(data, classes, missing_as=NA) {
  # where there is NA for a whole class of features, remove feature column
  if (is.na(missing_as)) {
    data[is.na(data)] = 0
  } else {
    data[which(data == missing_as, arr.ind=TRUE)] = 0
  }
  print(dim(data))
  print("Dropping features where at least one class is NA")
  uniq = unique(classes)
  subsets = list()
  for (i in uniq) {subsets[[i]] = data[grep(i, rownames(data)), ]}
  subsets = lapply(lapply(lapply(subsets, colSums), data.frame), t)
  subsets = do.call("rbind", subsets)
  rownames(subsets) = uniq

  if (is.na(missing_as)) {
    subsets[which(subsets == 0, arr.ind=TRUE)] = NA
  } else {
    subsets[which(subsets == 0, arr.ind=TRUE)] = missing_as
  }

  subsets = t(na.omit(t(subsets)))
  data = data[, c(colnames(subsets))]

  if (is.na(missing_as)) {
    data[which(data == 0, arr.ind=TRUE)] = NA
  } else {
    data[which(data == 0, arr.ind=TRUE)] = missing_as
  }
  return(data)
}

cor_imputed_unimputed = function(pca_withna, pca_impute, names) {
  # plots a heatmap of correlations: -> df, df, name
  x <- pca_withna
  y <- pca_impute
  z <- names
  print("Plotting correlation between unimputed and imputed components")
  print(ggplot(melt(cor(x$variates$X, y$variates$X)),
    aes(Var1, Var2, fill=value)) +
    ggtitle(paste(z, "correlation between imputed and unimputed data")) +
    geom_tile() +
    scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0,
      limit=c(-1,1), space="Lab", name="Pearson\nCorrelation") +
    theme_minimal()
  )
}

make_case_study_1 <- function() {
  infile_path_1 <- "../results/case_study_1/RData.RData"
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
  infile_path_1 <- "../results/case_study_1/RData.RData"
  unimputed_prot_path <- "../data/case_study_1/diablo_proteome.txt"
  unimputed_prot <- read.table(unimputed_prot_path, sep="\t", header=TRUE, row.names=1)
  unimputed_tran_path <- "../data/case_study_1/diablo_translatome.txt"
  unimputed_tran <- read.table(unimputed_tran_path, sep="\t", header=TRUE, row.names=1)
  outfile_dir <- "../results/manuscript_figures/"

  # case study 1
  load(infile_path_1)
  
  pdf(paste(outfile_dir, "case_1_extra_plots.pdf", sep=""))
  
  unimputed_prot[unimputed_prot == 0] <- NA
  unimputed_tran[unimputed_tran == 0] <- NA
  na_prop_prot <- show_na_prop(
    unimputed_prot, "Proteome"
  )
  na_prop_tran <- show_na_prop(
    unimputed_tran, "Translatome"
  )

  # need to drop col where all missing
  unimputed_prot <- remove_na_class(unimputed_prot, classes, missing_as=NA)
  unimputed_tran <- remove_na_class(unimputed_tran, classes, missing_as=NA)

  pca_unimputed_prot <- pca(unimputed_prot, ncomp=10)
  pca_unimputed_tran <- pca(unimputed_tran, ncomp=10)
  pca_imputed_prot <- pca(data$proteome, ncomp=10)
  pca_imputed_tran <- pca(data$translatome, ncomp=10)
  
  prot <- "Proteome"
  tran <- "Translatome"
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

make_case_study_1_multilevel <- function() {
  # case study 1 multilevel plots
  infile_path_1 <- "../results/case_study_1/RData.RData"
  outfile_dir <- "../results/manuscript_figures/"

  # case study 1
  load(infile_path_1)

  pdf(paste(outfile_dir, "case_1_multilevel.pdf", sep=""))

  data_prot <- pca(data_imp$proteome, ncomp=10, scale=TRUE, center=TRUE)
  plotIndiv(
    data_prot, ind.names=FALSE, group=classes, legend=TRUE,
    pch=pch, title=paste("Proteome (uncorrected)", "1/2"), comp=c(1,2),
    ellipse=FALSE, background=FALSE
  )
  plotIndiv(
    data_pca_multilevel$proteome, ind.names=FALSE, group=classes, legend=TRUE,
    pch=pch, title=paste("Proteome (corrected)", "1/2"), comp=c(1,2), 
    ellipse=FALSE, background=FALSE
  )

  data_tran <- pca(data_imp$translatome, ncomp=10, scale=TRUE, center=TRUE)
  plotIndiv(
    data_tran, ind.names=FALSE, group=classes, legend=TRUE,
    pch=pch, title=paste("Translatome (uncorrected)", "1/2"), comp=c(1,2),
    ellipse=FALSE, background=FALSE
  )
  plotIndiv(
    data_pca_multilevel$translatome, ind.names=FALSE, group=classes, legend=TRUE,
    pch=pch, title=paste("Translatome (corrected)", "1/2"), comp=c(1,2),
    ellipse=FALSE, background=FALSE
  )

  dev.off()

  # clear environment
  rm(list=ls())
}

make_case_study_1_diablo <- function () {
  infile_path_1 <- "../results/case_study_1/RData.RData"
  outfile_dir <- "../results/manuscript_figures/"

  # case study 1
  load(infile_path_1)
  pdf(paste(outfile_dir, "case_1_diablo.pdf", sep=""))
  circosPlot(diablo, cutoff=0.95, size.variables=0.01, size.labels=2)
  dev.off()
}

make_case_study_1_error <- function() {
  infile_path_1 <- "../results/case_study_1/RData.RData"
  outfile_dir <- "../results/manuscript_figures/"

  # case study 1
  load(infile_path_1)
  pdf(paste(outfile_dir, "case_1_error.pdf", sep=""))

  proteome_plsda <- data_plsda[,1]
  class(proteome_plsda) <- c("mixo_plsda", "mixo_pls", "DA")
  translatome_plsda <- data_plsda[,2]
  class(translatome_plsda) <- c("mixo_plsda", "mixo_pls", "DA")
  proteome_splsda <- data_splsda[,1]
  class(proteome_splsda) <- c("mixo_splsda", "mixo_spls", "DA")
  translatome_splsda <- data_splsda[,2]
  class(translatome_splsda) <- c("mixo_splsda", "mixo_spls", "DA")
  
  perf_plsda_prot <- perf(
    object=proteome_plsda, validation="loo", nrepeats=10, 
    auc=TRUE, cpus=6, progressBar=TRUE
  )
  perf_plsda_tran <- perf(
    object=translatome_plsda, validation="loo", nrepeats=10,
    auc=TRUE, cpus=6, progressBar=TRUE
  )
  perf_splsda_prot <- perf(
    object=proteome_splsda, validation="loo", nrepeats=10,
    auc=TRUE, cpus=6, progressBar=TRUE
  )
  perf_splsda_tran <- perf(
    object=translatome_splsda, validation="loo", nrepeats=10,
    auc=TRUE, cpus=6, progressBar=TRUE
  )
  perf_diablo <- perf(
    object=diablo, validation="loo", nrepeats=10,
    auc=TRUE, cpus=6, progressBar=TRUE
  )

  error_plsda <- c(
    apply(perf_plsda_prot$error.rate$overall, 1, min),
    apply(perf_plsda_tran$error.rate$overall, 1, min)
  )
  error_splsda <- c(
    apply(perf_splsda_prot$error.rate$overall, 1, min),
    apply(perf_splsda_tran$error.rate$overall, 1, min)
  )
  error_singleomics <- data.frame(c(
    as.vector(unlist(error_plsda)),
    as.vector(unlist(error_splsda))
  ))
  error_diablo <- min(
    perf_diablo$WeightedVote.error.rate$mahalanobis.dist["Overall.ER",]
  )

  # temporary variable for plotting only
  colnames(error_singleomics) <- "es"
  es <- error_singleomics
  es <- cbind(x=as.factor(seq(1)), es)
  es_plot <- ggplot(es, aes(x=x, y=es)) +
    geom_boxplot(color="orange", width=.3) +
    theme(axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
    scale_x_discrete(labels = NULL, breaks = NULL) + 
    theme_minimal() +
    ggtitle("Minimum error rates across all single omics data") +
    geom_hline(aes(yintercept=error_diablo), color="blue") +
    scale_y_continuous(
      breaks=sort(c(seq(min(es$es), max(es$es), length.out=4), error_diablo))
    ) +
    geom_text(
      aes(0,error_diablo,label="DIABLO minimum error rate", hjust=0, vjust=-1)
      ) +
    xlab("Case study 1") +
    ylab("Minimum error rates")
  print(es_plot)
  dev.off()

  # clear environment
  rm(list=ls())
}

make_case_study_2 <- function() {
  
  infile_path_2 <- "../results/case_study_2/RData.RData"
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
  infile_path_2 <- "../results/case_study_2/RData.RData"
  unimputed_path <- "../data/case_study_2/data_transcriptomics.tsv"
  unimputed <- read.table(unimputed_path, sep="\t", header=TRUE, row.names=1)
  outfile_dir <- "../results/manuscript_figures/"

  # case study 2
  load(infile_path_2)

  pdf(paste(outfile_dir, "case_2_extra_plots.pdf", sep=""))

  unimputed[unimputed == 0] <- NA
  na_prop <- show_na_prop(
    unimputed, "Transcriptome"
  )
  #unimputed <- remove_na_class(unimputed, classes, missing_as=NA)
  pca_imputed <- pca(data$transcriptome, ncomp=10)
  pca_unimputed <- pca(unimputed, ncomp=10)
  
  name <- "Transcriptome"
  cor_imputed_unimputed(pca_imputed, pca_unimputed, name)
  dev.off()

  # clear environment
  rm(list=ls())
}

make_case_study_2_diablo <- function () {
  infile_path_2 <- "../results/case_study_2/RData.RData"
  outfile_dir <- "../results/manuscript_figures/"

  # case study 1
  load(infile_path_2)
  pdf(paste(outfile_dir, "case_2_diablo.pdf", sep=""))
  circosPlot(diablo, cutoff=0.8, size.variables=0.01, size.labels=2)
  dev.off()
}

make_case_study_2_error <- function() {
  infile_path_2 <- "../results/case_study_2/RData.RData"
  outfile_dir <- "../results/manuscript_figures/"

  # case study 2
  load(infile_path_2)
  pdf(paste(outfile_dir, "case_2_error.pdf", sep=""))
  
  error_plsda <- c(
    apply(data_plsda$lipidome$perf_plsda$error.rate$overall, 1, min),
    apply(data_plsda$metabolome$perf_plsda$error.rate$overall, 1, min),
    apply(data_plsda$proteome$perf_plsda$error.rate$overall, 1, min),
    apply(data_plsda$transcriptome$perf_plsda$error.rate$overall, 1, min)
  )
  error_splsda <- c(
    apply(data_splsda$lipidome$perf_splsda$error.rate$overall, 1, min),
    apply(data_splsda$metabolome$perf_splsda$error.rate$overall, 1, min),
    apply(data_splsda$proteome$perf_splsda$error.rate$overall, 1, min),
    apply(data_splsda$transcriptome$perf_splsda$error.rate$overall, 1, min)
  )
  error_singleomics <- data.frame(c(
    as.vector(unlist(error_plsda)), 
    as.vector(unlist(error_splsda))
  ))
  error_diablo <- min(
    perf_diablo$WeightedVote.error.rate$centroids.dist["Overall.ER",]
  )
  # temporary variable for plotting only
  colnames(error_singleomics) <- "es"
  es <- error_singleomics
  es <- cbind(x=as.factor(seq(1)), es)
  es_plot <- ggplot(es, aes(x=x, y=es)) +
    geom_boxplot(color="orange", width=.3) +
    theme(axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    theme_minimal() +
    ggtitle("Minimum error rates across all single omics data") +
    geom_hline(aes(yintercept=error_diablo), color="blue") +
    scale_y_continuous(
      breaks=sort(c(seq(min(es$es), max(es$es), length.out=4), error_diablo))
    ) +
    geom_text(
      aes(0,error_diablo,label="DIABLO minimum error rate", hjust=0, vjust=-1)
      ) +
    xlab("Case study 2") + 
    ylab("Minimum error rates")
  print(es_plot)
  dev.off()
  
  # clear environment
  rm(list=ls())
}

main <- function() {
  # case study 1 was generated with an old version of the code
  # some data structure change but underlying info is the same
  make_case_study_1()
  make_case_study_1_extra()  
  make_case_study_1_multilevel()
  make_case_study_1_diablo()
  make_case_study_1_error()

  # case study 2 was generated with an up to date code version
  make_case_study_2()
  make_case_study_2_extra()
  make_case_study_2_diablo()
  make_case_study_2_error()
}

main()
