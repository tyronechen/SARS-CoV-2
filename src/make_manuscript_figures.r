#!/usr/bin/Rscript
# regenerate plots for this manuscript
library(mixOmics)

make_indiv <- function(data, classes, title, dist="centroids.dist", 
		       pch=NA, bg=FALSE) {
  # take a mixomics data object and make plots
  if (bg == TRUE) {
    bg = background.predict(data, comp.predicted=2, dist=dist)
  } else { bg <- NULL }
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

make_auroc <- function(data, title="", pch=NULL) {
  # take tuned mixomics object and plot roc curve
  data <- auroc(data, title=paste("ROC comp 1", title),
  		multilevel=pch)
  return(data)
}

make_all <- function(data, title, classes, outfile_path, 
		     pch=NULL, cpus=6, dist="centroids.dist", 
		     show_cols=FALSE, tuned=NA, bg=FALSE,
		     stability=FALSE, force_short=FALSE) {
  # aggregate all functions into master plotting function
  pdf(outfile_path)
  make_indiv(data=data, classes=classes, dist=dist,
             title=title, pch=pch, bg=bg)
  make_cim(data=data, classes=classes, show_cols=show_cols,
           title=title, pch=pch)
  make_loadings(data=data, title=title, force_short=force_short)
  make_error(data=data, title=title, stability=stability)
  if (!is.na(tuned)) {make_tuned(data=tuned, title=title)}
  make_auroc(data=data, title=title, pch=pch)
  dev.off()
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

main <- function() {
  infile_path_1 <- "../results/RData.RData"
  infile_path_2 <- ""
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

  # case study 2
  #load(infile_path_2)
  #load("./RData.RData")
}

main()
