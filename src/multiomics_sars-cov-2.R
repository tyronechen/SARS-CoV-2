#!/usr/bin/Rscript
library(argparser, quietly=TRUE)
library(igraph)
library(mixOmics)

parse_data = function(infile_path, offset=0, missing_as=NA, rmna=TRUE) {
  # load in omics data into a diablo-compatible format: infile_path -> dataframe
  print("Parsing file:")
  print(infile_path)
  data = read.table(infile_path, sep="\t", header=TRUE, row.names=1) + offset
  if (is.na(missing_as)) {
    data[which(data == 0, arr.ind=TRUE)] = NA
  } else {
    data[which(data == 0, arr.ind=TRUE)] = missing_as
  }
  if (rmna == TRUE) {
    print("Dropping features where all values are NA")
    data = data[, unlist(lapply(data, function(x) !all(is.na(x))))]
  }
  return(data)
}

parse_mappings = function(infile_path) {
  # load in id:name data into a mapping table, format: infile_path -> dataframe
  print("Parsing mappings...")
  print(infile_path)
  data = read.table(infile_path, sep="\t", header=TRUE, row.names=2)
  data["X"] = NULL
  return(data)
}

remap_data = function(data, mapping) {
  # map feature names to feature id, format: dataframe, dataframe -> dataframe
  # TODO: doesnt handle duplicates in row ids (some feature id map to same name)
  print("Remapping data with map files (duplicate row names will be renamed!)")
  data = merge(t(data), mapping, by=0, all.x=TRUE)
  row.names(data) = make.names(data$val, unique=TRUE)
  data["Row.names"] = NULL
  data["val"] = NULL
  return(t(data))
}

show_na_prop = function(data_na, name) {
  # first, sum the number of missing values per variable
  sum_na_per_var = apply(data_na, 2, function(x) {sum(is.na(x))})
  # show proportion of NA values across all samples (y) for a variable (x)
  plot(sum_na_per_var/nrow(data_na), type='h', xlab='variable index',
    ylab='NA proportion (across all samples for a variable)',
    main=paste(name, 'NA rate per variable on unfiltered data'))
}

remove_na_prop = function(data, class, pch=NA, na_prop=0.3) {
  mapply(function(x) remove_na_prop_(x, class, pch, na_prop), data)
}

remove_na_prop_ = function(data, class, pch=NA, na_prop=0.3) {
  # we want to compare dimensions later to see how many features were dropped
  data_na = data

  # first, sum the number of missing values per variable
  sum_na_per_var = apply(data_na, 2, function(x){sum(is.na(x))})

  # these variables could be removed
  remove_var = which(sum_na_per_var/ncol(data_na) >= na_prop)
  if (length(remove_var) > 0) {
    print("Removing this many columns from data:")
    print(length(remove_var))
    data_na = data_na[, -c(remove_var)]
  } else {
    print("No columns removed from data")
  }

  print("Data dimensions (original):")
  print(dim(data))
  print("Data dimensions (filtered):")
  print(dim(data_na))

  print("Proportion of missing values (original):")
  print(sum(is.na(data)) / (nrow(data) * ncol(data)))
  print("Proportion of missing values (filtered):")
  print(sum(is.na(data_na)) / (nrow(data_na) * ncol(data_na)))

  return(data_na)
  q()
  # do imputation
  nipals_X = nipals(data_na, reconst = TRUE, ncomp = 3)$rec

  # replace NA values only with imputed values
  id_na = is.na(data_na)
  nipals_X[!id_na] = data_na[!id_na]

  # compare pcas before and after imputation
  if (!is.na(pch)) {
    pca_with_na <- pca(data_na, ncomp=3, center=TRUE, scale=TRUE, pch=pch)
    pca_no_na <- pca(nipals_X, ncomp=3, center=TRUE, scale=TRUE, pch=pch)
  } else {
    pca_with_na <- pca(data_na, ncomp=3, center=TRUE, scale=TRUE, pch=pch)
    pca_no_na <- pca(nipals_X, ncomp=3, center=TRUE, scale=TRUE, pch=pch)
  }

  print(pca_with_na$cum.var)
  print(pca_no_na$cum.var)
  plotIndiv(pca_with_na, group = class, legend = TRUE)
  plotIndiv(pca_no_na, group = class, legend = TRUE)

  plot(cor(pca_with_na$variates$X, pca_no_na$variates$X))
}

remove_na_class = function(data, classes, missing_as=NA) {
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

remove_novar = function(data) {
  # samples with zero variance are meaningless for PCA: dataframe -> dataframe
  # print("Dimensions before removing invariant columns:")
  # print(dim(data))
  data = data[, which(apply(data, 2, var) != 0)]
  # print("Dimensions after removing invariant columns:")
  # print(dim(data))
  return(data)
}

parse_classes = function(infile_path) {
  # load in class data for diablo: infile_path -> vector (of strings)
  data = read.table(infile_path, sep="\t", header=TRUE, row.names=1)
  return(unlist(as.vector(t(data))))
}

create_design = function(data, link=0.1) {
  # create design matrix from data: dataframe, link -> matrix (design)
  design = matrix(link, ncol = length(data), nrow = length(data),
                  dimnames = list(names(data), names(data)))
  diag(design) = 0
  return(design)
}

zero_to_na = function(data) {
  data[which(data == 0, arr.ind=TRUE)] = NA
  return(data)
}

count_missing = function(data) {
  # show missing values in data: dataframe -> list(dataframe, double)
  # return 0 values to NA
  ids_na = is.na(data)
  pct_na = sum(is.na(data)) / (nrow(data) * ncol(data))
  print("Percentage of missing values in data:")
  print(pct_na)
  return(list(ids_na=ids_na, pct_na=pct_na))
}

impute_missing_ = function(data, ncomp=10, block_name="", outdir="./") {
  # impute missing values with nipals: dataframe (with NA) -> dataframe
  print("Number of components for imputation:")
  print(ncomp)
  # data <- data[!is.na(apply(data, 1, function(x) var(x, na.rm = TRUE))),
  #              !is.na(apply(data, 2, function(x) var(x, na.rm = TRUE)))]
  nipals_tune = nipals(data, reconst=TRUE, ncomp=ncomp)
  barplot(nipals_tune$eig, main=paste(block_name, "Screeplot (nipals imputed)"),
    xlab="Number of components", ylab="Explained variance"
  )
  outfile_path = paste(outdir, "/", "data_", block_name, "_imputed.tsv", sep="")
  write.table(
    as.data.frame(nipals_tune$rec), file=outfile_path, quote=FALSE, sep="\t"
  )
  return(nipals_tune$rec)
}

impute_missing = function(data, ncomps, outdir) {
  mapply(function(x, y, z)
    impute_missing_(x, y, z, outdir), data, ncomps, names(data), SIMPLIFY=FALSE
  )
}

replace_missing_ = function(data, imputed) {
  mask = is.na(data)
  imputed[!mask] = data[!mask]
  return(imputed)
}

replace_missing = function(data, imputed) {
  mapply(function(x, y) replace_missing_(x, y), data, imputed)
}

plot_pca_single = function(data, classes, pch=NA, title="", ncomp=0, show=FALSE) {
  # do pca on individual classes: dataframe, vector, vector -> outfile_path.pdf
  names = names(data)

  print("Removing 0 variance columns from data...")
  data = lapply(data, remove_novar)

  if (ncomp == 0) {ncomp = dim(classes)[1]}

  data_pca = lapply(data, pca, ncomp=ncomp, center=TRUE, scale=TRUE)

  if (show == TRUE) {
    print("Showing PCA component contribution...")
    print(data_pca)
  }

  print("Plotting PCA component contribution...")
  mapply(function(x, y) plot(x, main=paste(y, "Screeplot")), data_pca, names)

  if (!is.na(pch)) {
    print("Plotting PCA by groups...")
    mapply(function(x, y) plotIndiv(x, comp=c(1,2), ind.names=FALSE,
      group=classes, legend=TRUE, ncomp=ncomp,
      title=paste(title, y, "PCA 1/2"), pch=pch), data_pca, names)
    mapply(function(x, y) plotIndiv(x, comp=c(1,3), ind.names=FALSE,
      group=classes, legend=TRUE, ncomp=ncomp,
      title=paste(title, y, "PCA 1/3"), pch=pch), data_pca, names)
    mapply(function(x, y) plotIndiv(x, comp=c(2,3), ind.names=FALSE,
      group=classes, legend=TRUE, ncomp=ncomp,
      title=paste(title, y, "PCA 2/3"), pch=pch), data_pca, names)
  } else {
    print("Plotting PCA by groups...")
    mapply(function(x, y) plotIndiv(x, comp=c(1,2), ind.names=TRUE,
      group=classes, legend=TRUE, ncomp=ncomp,
      title=paste(title, y, "PCA 1/2")), data_pca, names)
    mapply(function(x, y) plotIndiv(x, comp=c(1,3), ind.names=TRUE,
      group=classes, legend=TRUE, ncomp=ncomp,
      title=paste(title, y, "PCA 1/3"), pch=pch), data_pca, names)
    mapply(function(x, y) plotIndiv(x, comp=c(2,3), ind.names=TRUE,
      group=classes, legend=TRUE, ncomp=ncomp,
      title=paste(title, y, "PCA 2/3"), pch=pch), data_pca, names)
  }
  return(data_pca)
}

plot_additional = function(data, data_pca, names) {
  # correlation circle and biplots: dataframe, list, vector -> outfile_path.pdf
  print("Plotting correlation circle plots...")
  mapply(function(x, y) plotVar(x, comp=c(1, 2), title=paste(y, "PCA 1/2"),
    var.names=FALSE),
    data_pca, names)

  print("Plotting biplots...")
  mapply(function(x, y, z) biplot(y, cex=0.7, xlabs=paste(classes, 1:nrow(x)),
    main=paste(z, "Biplot")), data, data_pca, names)
}

plot_pca_multilevel = function(data, classes, pch, title="", ncomp=0, show=FALSE) {
  names = names(data)

  print("Removing 0 variance columns from data...")
  data = lapply(data, remove_novar)

  if (ncomp == 0) {ncomp = dim(classes)[1]}

  data_pca = lapply(data,pca,ncomp=ncomp,center=TRUE,scale=TRUE,multilevel=pch)
  if (show == TRUE) {
    print("Showing PCA multilevel component contribution...")
    print(data_pca)
  }

  print("Plotting PCA multilevel component contribution...")
  mapply(function(x, y) plot(x, main=paste(y, "Screeplot multilevel")),
    data_pca, names)

  print("Plotting PCA multilevel...")
  mapply(function(x, y) plotIndiv(x, comp=c(1,2), ind.names=FALSE,
    group=classes, legend=TRUE, ncomp=ncomp,
    title=paste(title, y, "PCA M 1/2"), pch=pch), data_pca, names)
  mapply(function(x, y) plotIndiv(x, comp=c(1,3), ind.names=FALSE,
    group=classes, legend=TRUE, ncomp=ncomp,
    title=paste(title, y, "PCA M 1/3"), pch=pch), data_pca, names)
  mapply(function(x, y) plotIndiv(x, comp=c(2,3), ind.names=FALSE,
    group=classes, legend=TRUE, ncomp=ncomp,
    title=paste(title, y, "PCA M 2/3"), pch=pch), data_pca, names)
  return(data_pca)
}

classify_plsda = function(data, classes, pch=NA, title="", ncomp=0,
  contrib="max", outdir="./", mappings=NULL, dist="centroids.dist", bg=TRUE) {
  mapply(function(x, y) classify_plsda_(
    x, classes, pch, y, ncomp, contrib, outdir, bg), data, title)
}

classify_plsda_ = function(data, classes, pch=NA, title="", ncomp=0,
  contrib="max", outdir="./", mappings=NULL, dist="centroids.dist", bg=TRUE) {
  # discriminate samples: list, vector, bool, integer -> list
  # single or multilevel PLS-DA
  if (!is.na(pch)) {
    print("Plotting multi level partial least squares discriminant analysis")
    pch = c(as.factor(pch))
    title_plt = paste(title, "PLSDA multi")
    data_plsda = plsda(data, Y=classes, multilevel=c(as.factor(pch)), ncomp=ncomp)
    if (!is.na(bg)) {
      bg = background.predict(data_plsda, comp.predicted=2, dist=dist)
      plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
        pch=pch, title=paste(title_plt, "1/2"), comp=c(1,2), ellipse=TRUE,
        background=bg
      )
    }
    plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
      pch=pch, title=paste(title_plt, "1/2"), comp=c(1,2), ellipse=TRUE,
    )
    plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
      pch=pch, title=paste(title_plt, "1/3"), comp=c(1,3), ellipse=TRUE
    )
    plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
      pch=pch, title=paste(title_plt, "2/3"), comp=c(2,3), ellipse=TRUE
    )
  } else {
    print("Plotting multi level partial least squares discriminant analysis")
    title_plt = paste(title, "PLSDA single")
    data_plsda = plsda(data, Y=classes, ncomp=ncomp)
    if (!is.na(bg)) {
      bg = background.predict(data_plsda, comp.predicted=2, dist=dist)
      plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
        title=paste(title_plt, "1/2"), comp=c(1,2), ellipse=TRUE,
        background=bg
      )
    }
    plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
      title=paste(title_plt, "1/2"), comp=c(1,2), ellipse=TRUE,
    )
    plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
      title=paste(title_plt, "1/3"), comp=c(1,3), ellipse=TRUE
    )
    plotIndiv(data_plsda, ind.names=FALSE, group=classes, legend=TRUE,
      title=paste(title_plt, "2/3"), comp=c(2,3), ellipse=TRUE
    )
  }

  print("Getting performance metrics")
  print("Plotting error rates...")
  metrics = perf(data_plsda, validation="loo", progressBar=TRUE, auc=TRUE)
  plot(metrics, main="Error rate PLSDA", col=color.mixo(5:7), sd=TRUE)

  # supporess roc printing automatically to stdout
  sink("/dev/null")
  roc = mapply(function(x) auroc(data_plsda, roc.comp=x), seq(ncomp))
  sink()

  # consider using plotvar for some datasets
  if (ncomp > 1) {
    print("Plotting arrow plot...")
    plotArrow(data_plsda, ind.names=FALSE, legend=TRUE, title="PLSDA")
  }
  # plotVar(data_plsda, legend=TRUE)

  print("Getting loadings and plotting clustered image maps")

  # setup colour map for clustered image plots
  colours_class = color.mixo(1:length(unique(classes)))[as.numeric(as.factor(classes))]

  if (!is.na(pch)) {
    colours_pch = color.mixo(1:length(unique(pch)))[as.numeric(as.factor(pch))]
    colours_cim = cbind(colours_class, colours_pch)
  } else {colours_cim = data.frame(colours_class)}

  # hide column (sample) names by default (will run off page otherwise)
  if (max(unlist(lapply(data_plsda$names$colnames$X, nchar))) > 8) {
    show_cols = FALSE
  } else {
    show_cols = TRUE
  }

  cim(data_plsda, title="PLSDA", row.sideColors=colours_cim,
    legend=list(title="Status")
  )
  for (comp in seq(ncomp)) {
    cim(data_plsda, comp=comp, title=paste("PLSDA Component", comp),
      row.sideColors=colours_cim, legend=list(title="Status")
    )
    plotLoadings(data_plsda, contrib="max", comp=comp, max.name.length=16,
      method='median', ndisplay=20, name.var=colnames(data), size.name=0.6,
      size.legend=0.6, title=paste(title, comp, "PLSDA max loadings"))
    plotLoadings(data_plsda, contrib="min", comp=comp, max.name.length=16,
      method='median', ndisplay=20, name.var=colnames(data), size.name=0.6,
      size.legend=0.6, title=paste(title, comp, "PLSDA min loadings"))
    loading_max = plotLoadings(data_plsda, contrib="max", comp=comp,
      method='median', ndisplay=NULL, name.var=colnames(data), plot=FALSE)
    loading_min = plotLoadings(data_plsda, contrib="min", comp=comp,
      method='median', ndisplay=NULL, name.var=colnames(data), plot=FALSE)
    title = gsub(" ", "_", title)
    path_max = paste(outdir, "/", title, "_", comp, "_PLSDA_max.txt", sep="")
    path_min = paste(outdir, "/", title, "_", comp, "_PLSDA_min.txt", sep="")
    print("Writing PLSDA loadings to:")
    print(path_max)
    print(path_min)
    write.table(as.data.frame(loading_max), file=path_max, quote=FALSE, sep="\t")
    write.table(as.data.frame(loading_min), file=path_min, quote=FALSE, sep="\t")
  }
  print(metrics)
  return(list(data_plsda=data_plsda, perf_plsda=metrics))
}

tune_splsda = function(data, classes, names, multilevel=NULL, ncomp=3, nrepeat=10,
  logratio="none", test_keepX=c(5, 50, 100), validation="loo", folds=10,
  dist="centroids.dist", cpus=2, progressBar=TRUE) {
    mapply(function(x, y) tune_splsda_(x, classes, names, multilevel, ncomp,
      nrepeat, logratio, test_keepX, validation, folds, dist, cpus, progressBar),
      data, names, SIMPLIFY=FALSE)
  }

tune_splsda_ = function(data, classes, names, multilevel=NULL, ncomp=0, nrepeat=10,
  logratio="none", test_keepX=c(5, 50, 100), validation="loo", folds=10,
  dist="centroids.dist", cpus=2, progressBar=TRUE) {
  if (ncomp == 0) {ncomp = (length(test_keepX))}
  # tune splsda components
  tuned = tune.splsda(data, Y=classes, multilevel=multilevel, ncomp=ncomp,
    nrepeat=nrepeat, logratio=logratio, test.keepX=test_keepX,
    validation=validation, folds=folds, dist=dist, cpus=cpus,
    progressBar=progressBar
  )
  print(plot(tuned, main=names))
  return(tuned)
}

classify_splsda = function(data, classes, pch=NA, title="", ncomp=NULL,
  keepX=NULL, contrib="max", outdir="./", mappings=NULL, dist="centroids.dist",
  bg=TRUE) {
  mapply(function(x, y, c, k) classify_splsda_(
    x, classes, pch, y, c, k, contrib, outdir
  ), data, title, ncomp, keepX)
}

classify_splsda_ = function(data, classes, pch=NA, title="", ncomp=NULL,
  keepX=NULL, contrib="max", outdir="./", mappings=NULL, dist="centroids.dist",
  bg=TRUE) {
  # discriminate samples: list, vector, bool, integer, vector -> list
  # single or multilevel sPLS-DA
  if (is.null(keepX)) {
    print("The number of variables selected on each component is not selected!")
    q(status=1)
  }
  if (is.null(ncomp)) {
    print("Invalid number of components, inferring from keepX")
    ncomp = length(keepX)
    print(ncomp)
  }
  print("splsda components:")
  print(ncomp)
  print("number of variables on each component:")
  print(keepX)
  classes = as.factor(classes)

  if (!is.na(pch)) {
    data_splsda = splsda(data,Y=classes,multilevel=pch,ncomp=ncomp,keepX=keepX)
  } else {
    data_splsda = splsda(data, Y=classes, ncomp=ncomp, keepX=keepX)
  }

  if (!is.na(bg)) {
    bg = background.predict(data_splsda, comp.predicted=2, dist=dist)
    plotIndiv(data_splsda, ind.names=FALSE, group=classes, legend=TRUE,
      pch=pch, title=paste(title, "sPLSDA multi 1/2"), comp=c(1,2),
      ellipse=TRUE, background=bg
    )
  } else {
    plotIndiv(data_splsda, ind.names=FALSE, group=classes, legend=TRUE,
      pch=pch, title=paste(title, "sPLSDA multi 1/2"), comp=c(1,2),
      ellipse=TRUE, background=NULL
    )
  }

  plotIndiv(data_splsda, ind.names=FALSE, group=classes, legend=TRUE,
    pch=pch, title=paste(title, "sPLSDA multi 1/2"), comp=c(1,2), ellipse=TRUE
  )
  plotIndiv(data_splsda, ind.names=FALSE, group=classes, legend=TRUE,
    pch=pch, title=paste(title, "sPLSDA multi 1/3"), comp=c(1,3), ellipse=TRUE
  )
  plotIndiv(data_splsda, ind.names=FALSE, group=classes, legend=TRUE,
    pch=pch, title=paste(title, "sPLSDA multi 2/3"), comp=c(2,3), ellipse=TRUE
  )

  print("Getting performance metrics")
  print("Plotting error rates...")
  metrics = perf(data_splsda, validation="loo", progressBar=TRUE, auc=TRUE)
  print(metrics$error.rate)
  plot(metrics, main="Error rate sPLSDA", col=color.mixo(5:7), sd=TRUE)
  print("Plotting stability of sPLSDA...")
  plot(metrics$features$stable[[1]], type="h", main="Comp 1", las=2,
    ylab="Stability", xlab="Features", xaxt='n'
  )
  plot(metrics$features$stable[[2]], type="h", main="Comp 2", las=2,
    ylab="Stability", xlab="Features", xaxt='n'
  )
  plot(metrics$features$stable[[3]], type="h", main="Comp 3", las=2,
    ylab="Stability", xlab="Features", xaxt='n'
  )
  sink("/dev/null")
  roc = mapply(function(x) auroc(data_splsda, roc.comp=x), seq(ncomp))
  sink()

  if (ncomp > 1) {
    print("Plotting arrow plot...")
    plotArrow(data_splsda, ind.names=FALSE, legend=TRUE, title="sPLSDA")
  }
  # plotVar(data_plsda, legend=TRUE)

  print("Getting loadings and plotting clustered image maps")

  # setup colour map for clustered image plots
  colours_class = color.mixo(1:length(unique(classes)))[as.numeric(as.factor(classes))]

  if (!is.na(pch)) {
    colours_pch = color.mixo(1:length(unique(pch)))[as.numeric(as.factor(pch))]
    colours_cim = cbind(colours_class, colours_pch)
  } else {colours_cim = data.frame(colours_class)}

  # hide column (sample) names by default (will run off page otherwise)
  if (max(unlist(lapply(data_splsda$names$colnames$X, nchar))) > 8) {
    show_cols = FALSE
  } else {
    show_cols = TRUE
  }
  cim(data_splsda, title="sPLSDA", row.sideColors=colours_cim,
    legend=list(title="Status"), col.names=show_cols
  )

  short = make.names(sapply(colnames(data), strtrim, 6, USE.NAMES=FALSE), unique=TRUE)
  for (comp in seq(ncomp)) {
    cim(data_splsda, comp=comp, title=paste("sPLSDA Component", comp),
      row.sideColors=colours_cim, legend=list(title="Status")
    )
    plotLoadings(data_splsda, contrib="max", comp=comp, max.name.length=8,
      method='median', ndisplay=20, name.var=short, size.name=0.6,
      size.legend=0.6, title=paste(title, comp, "sPLSDA max loadings"))
    plotLoadings(data_splsda, contrib="min", comp=comp, max.name.length=8,
      method='median', ndisplay=20, name.var=short, size.name=0.6,
      size.legend=0.6, title=paste(title, comp, "sPLSDA min loadings"))
    loading_max = plotLoadings(data_splsda, contrib="max", comp=comp,
      method='median', ndisplay=NULL, name.var=colnames(data), plot=FALSE)
    loading_min = plotLoadings(data_splsda, contrib="min", comp=comp,
      method='median', ndisplay=NULL, name.var=colnames(data), plot=FALSE)
    title = gsub(" ", "_", title)
    path_max = paste(outdir, "/", title, "_", comp, "_sPLSDA_max.txt", sep="")
    path_min = paste(outdir, "/", title, "_", comp, "_sPLSDA_min.txt", sep="")
    print("Writing sPLSDA loadings to:")
    print(path_max)
    print(path_min)
    write.table(as.data.frame(loading_max), file=path_max, quote=FALSE, sep="\t")
    write.table(as.data.frame(loading_min), file=path_min, quote=FALSE, sep="\t")
  }
  return(list(data_splsda=data_splsda, perf_splsda=metrics))
}

plot_plsda = function(data, classes, pch, title="", ncomp=0) {
  names = names(data)
  print("Plotting PLSDA component contribution...")
  mapply(function(x, y) plot(x, main=paste(y, "Screeplot multilevel")),
    data, names)

  print("Plotting plsda...")
  mapply(function(x, y) plotIndiv(x, comp=c(1,2), ind.names=TRUE,
    group=classes, legend=TRUE, ncomp=ncomp,
    title=paste(title, y, "PLSDA 1/2"), pch=pch), data_pca, names)
  mapply(function(x, y) plotIndiv(x, comp=c(1,3), ind.names=TRUE,
    group=classes, legend=TRUE, ncomp=ncomp,
    title=paste(title, y, "PLSDA 1/3"), pch=pch), data_pca, names)
  mapply(function(x, y) plotIndiv(x, comp=c(2,3), ind.names=TRUE,
    group=classes, legend=TRUE, ncomp=ncomp,
    title=paste(title, y, "PLSDA 2/3"), pch=pch), data_pca, names)
  return(data_pca)
}

plot_loadings = function(data, contrib="max", ncomp=2, method="median",
  ndisplay=50, title="Loadings") {
  mapply(function(x, y) plot_loadings_(x, contrib, ncomp, method, ndisplay,
    y), data, title)
}

plot_loadings_ = function(data, contrib="max", ncomp=2, method="median",
  ndisplay=50, title="Loadings") {
  name_var = names(data)
  loadings=plotLoadings(data, contrib, ncomp, method, ndisplay, name_var, title)
  print(loadings)
  path = paste(title, "txt", sep=".")
  print(path)
  write.table(as.data.frame(loadings), file=path, quote=FALSE, sep="\t")
  return(loadings)
}

cor_imputed_unimputed = function(pca_withna, pca_impute, names) {
  # plots a heatmap of correlations: -> list of df, list of df, vector of names
  print("Plotting correlation between unimputed and imputed components")
  mapply(function(x, y, z) print(ggplot(melt(cor(x$variates$X, y$variates$X)),
    aes(Var1, Var2, fill=value)) +
    ggtitle(paste(z, "Correlation between imputed and unimputed data")) +
    geom_tile() +
    scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0,
      limit=c(-1,1), space="Lab", name="Pearson\nCorrelation") +
    theme_minimal()),
  pca_withna, pca_impute, names)
}

tune_diablo_ncomp = function(data, classes, design, ncomp=0) {
  # First, we fit a DIABLO model without variable selection to assess the global
  # performance and choose the number of components for the final DIABLO model.
  # The function perf is run with 10-fold cross validation repeated 10 times.
  print("Finding optimal number of components for DIABLO...")
  if (ncomp == 0) {ncomp = length(unique(classes))}
  sgccda_res = block.splsda(X=data, Y=classes, ncomp=ncomp, design=design)

  # this code takes a couple of min to run
  perf_diablo = perf(sgccda_res, validation = 'loo', folds = 10, nrepeat = 10)

  # print(perf.diablo)  # lists the different outputs
  plot(perf_diablo, main="DIABLO optimal components")
  # perf_diablo$choice.ncomp$WeightedVote
  print(perf_diablo$choice.ncomp)

  print("Plotting stability of DIABLO components...")
  # return(perf_diablo)
  for (i in perf_diablo$features$stable$nrep1) {
    for (j in names(perf_diablo$features$stable$nrep1)) {
      if (any(is.infinite(i$comp1)) == FALSE) {
        plot(i$comp1, type="h", las=2, ylab="Stability", xlab="Features",
          main=paste(j, "Comp 1"), xaxt='n'
        )
        print(paste(j, "Comp 1"))
      }
      if (any(is.infinite(i$comp2)) == FALSE) {
        plot(i$comp2, type="h", las=2, ylab="Stability", xlab="Features",
          main=paste(j, "Comp 2"), xaxt='n'
        )
        print(paste(j, "Comp 2"))
      }
      if (any(is.infinite(i$comp3)) == FALSE) {
        plot(i$comp1, type="h", las=2, ylab="Stability", xlab="Features",
          main=paste(j, "Comp 3"), xaxt='n'
        )
        print(paste(j, "Comp 3"))
      }
    }
  }
  # sink("/dev/null")

  return(perf_diablo)
}

tune_diablo_keepx = function(data, classes, ncomp, design,
  test_keepX=c(5,50,100), cpus=2, dist="centroids.dist", progressBar=TRUE) {
  # This tuning function should be used to tune the keepX parameters in the
  #   block.splsda function.
  # We choose the optimal number of variables to select in each data set using
  # the tune function, for a grid of keepX values. Note that the function has
  # been set to favor the small-ish signature while allowing to obtain a
  # sufficient number of variables for downstream validation / interpretation.
  # See ?tune.block.splsda.
  print("Tuning keepX parameter...")
  test_keepX = mapply(function(name, dims) list(name=dims), names(data),
    rep(list(test_keepX))
  )

  tune_data = tune.block.splsda(
      X=data, Y=classes, ncomp=ncomp, test.keepX=test_keepX, design=design,
      validation='loo', folds=10, nrepeat=1, cpus=cpus, dist=dist,
      progressBar=progressBar)
  list_keepX = tune_data$choice.keepX
  return(list_keepX)
}

force_unique_blocks = function(data) {
  # in diablo, features across blocks must be unique: list of df -> list of df
  print("Appending suffix to individual block names (diablo requires unique!):")
  names = names(data)
  # print(names)
  colnames_new = mapply(
    function(x, y) paste(x, y, sep="_"), lapply(data, colnames), names(data)
  )
  reassign_colnames_ = function(data, colnames_new) {
    colnames(data) = colnames_new
    return(data)
  }
  data = mapply(reassign_colnames_, data, colnames_new)
  names(data) = names
  return(data)
}

run_diablo = function(data, classes, ncomp, design, keepx=NULL) {
  # this is the actual part where diablo is run
  print("Running DIABLO...")
  block.splsda(X=data, Y=classes, ncomp=ncomp, keepX=keepx, design=design)
}

plot_diablo = function(data, ncomp=0, outdir="./", data_names=NA, keepvar="") {
  # plot the diablo data with a series of diagnostic plots

  # need to make a function to squeeze sample names automatically and remap
  trim_names_ = function(data, trim=16) {
    all_names = data
    long_names = which(sapply(data, nchar, USE.NAMES=FALSE) > trim)
    if (length(long_names) == 0) {return()}
    original_names = data[long_names]
    for (i in long_names) { data[i] = as.character(i) }
    # later map these back
    maptable = data.frame(from=long_names, to=original_names)
    return(list(data=data, all_names=all_names, maptable=maptable))
  }
  trimmed_names = lapply(data$names$colnames, trim_names_)
  block_to_trim = names(trimmed_names[lapply(trimmed_names, length) > 0])

  # replace names in all associated columns for visualisation only
  replace_names_ = function(data, trim=16) {
    all_names = head(data$names$colnames, n=-1)
    split_ = function(block, all_names) {
      all_names = gsub("__FEATUREID", "", all_names[[block]])
      all_names = gsub(paste("_", block, sep=""), "", all_names)
      return(all_names)
    }
    splitted = mapply(
      function(x) split_(x, all_names), names(all_names), SIMPLIFY=FALSE
    )

    truncate_ = function(names, trim) {
      ifelse(nchar(names) > trim, paste0(strtrim(names, trim), ''), names)
    }
    truncated = lapply(splitted, truncate_, trim)

    # make a copy, dont want to overwrite
    data_vis = data

    # map truncated values to all locations
    for (i in names(all_names)) {
      row.names(data_vis$loadings[[i]]) = make.unique(sapply(truncated[[i]], toString), sep="__")
      data_vis$names$colnames[[i]] = make.unique(sapply(truncated[[i]], toString), sep="__")
      colnames(data_vis$X[[i]]) = make.unique(sapply(truncated[[i]], toString), sep="__")
    }
    return(list(data_vis=data_vis, truncated=truncated))
  }
  data_vis_names = replace_names_(data, trim=16)
  data_vis = data_vis_names$data_vis
  truncated = data_vis_names$truncated
  print("Plotting correlation between components...")
  # roc = mapply(function(x) auroc(data_plsda, roc.comp=x), seq(ncomp))
  mapply(function(x) plotDiablo(data, ncomp=x), seq(ncomp))
  # plotDiablo(data, ncomp = 1)
  if (ncomp > 1) {
    print("Plotting individual samples into space spanned by block components...")
    plotIndiv(data_vis, ind.names=FALSE, legend=TRUE, title='DIABLO', ellipse=TRUE)
    print("Plotting arrow plot...")
    plotArrow(data_vis, ind.names=FALSE, legend=TRUE, title='DIABLO')
  }
  print("Plotting correlation circle plot...")
  plotVar(data_vis, style='graphics', legend=TRUE, comp=c(1,2),
    title="DIABLO 1/2", var.names=FALSE
  )
  if (ncomp > 2) {
    plotVar(data_vis, style='graphics', legend=TRUE, comp=c(1,3),
      title="DIABLO 1/3", var.names=FALSE
    )
    plotVar(data_vis, style='graphics', legend=TRUE, comp=c(2,3),
      title="DIABLO 2/3", var.names=FALSE
    )
  }
  print("Plotting circos from similarity matrix...")
  # cant remove feature labels, need to make label size 0.001 or lower
  corr_diablo = circosPlot(
    data, cutoff=0.95, line=TRUE, size.legend=0.5, size.variables=0.001,
    var.names=truncated
  )
  corr_out = file=paste(outdir,"/DIABLO_var_",keepvar,"_correlations.txt",sep="")
  write.table(corr_diablo, file=corr_out, sep="\t", quote=FALSE)
  print("Plotting relevance network from similarity matrix...")
  cyto = network(
    data, blocks=c(1,2), color.node=c('darkorchid','lightgreen'), cutoff=0.4
  )
  cyto_out = paste(outdir, "/DIABLO_var_", keepvar, "_network.graphml", sep="")
  write.graph(cyto$gR, cyto_out, format="graphml")
  print("Plotting overall heatmap...")

  # hide column (sample) names by default (will run off page otherwise)
  if (max(unlist(lapply(data$names$colnames$X, nchar))) > 8) {
    show_cols = FALSE
  } else {
    show_cols = TRUE
  }

  cimDiablo(data, size.legend=0.5, col.names=show_cols)

  block_to_trim = names(trimmed_names[lapply(trimmed_names, length) > 0])

  print("Plotting loading weight of selected variables on each component...")
  for (comp in seq(ncomp)) {
    for (i in block_to_trim) {
      data$names$colnames[[i]] = trimmed_names[[i]][["data"]]
    }
    cimDiablo(data, comp=comp, size.legend=0.5)
    plotLoadings(data, contrib="max", comp=comp, max.name.length=8,
      method='median', ndisplay=20, name.var=colnames(data), size.name=0.6,
      size.legend=0.6, title=paste(comp, "DIABLO max loadings"))
    plotLoadings(data, contrib="min", comp=comp, max.name.length=8,
      method='median', ndisplay=20, name.var=colnames(data), size.name=0.6,
      size.legend=0.6, title=paste(comp, "DIABLO min loadings"))
    for (i in block_to_trim) {
      data$names$colnames[[i]] = trimmed_names[[i]][["all_names"]]
    }

    for (i in data_names) {
      for (i in block_to_trim) {
        data$names$colnames[[i]] = trimmed_names[[i]][["data"]]
      }
      plotLoadings(data, contrib="max", comp=comp, block=i, max.name.length=8,
        method='median', ndisplay=20, name.var=colnames(data), plot=TRUE,
        title=paste(comp, i, "DIABLO max loadings"), size.name=0.6
      )
      plotLoadings(data, contrib="min", comp=comp, block=i, max.name.length=8,
        method='median', ndisplay=20, name.var=colnames(data), plot=TRUE,
        title=paste(comp, i, "DIABLO min loadings"), size.name=0.6
      )
      for (i in block_to_trim) {
        data$names$colnames[[i]] = trimmed_names[[i]][["all_names"]]
      }

      loading_max = plotLoadings(data, contrib="max", comp=comp, block=i,
        method='median', ndisplay=NULL, name.var=colnames(data), plot=FALSE)
      loading_min = plotLoadings(data, contrib="min", comp=comp, block=i,
        method='median', ndisplay=NULL, name.var=colnames(data), plot=FALSE)
      # title = gsub(" ", "_", title)
      path_max = paste(
        outdir, "/", i, "_", comp, "_DIABLO_var_", keepvar, "_max.txt", sep=""
      )
      path_min = paste(
        outdir, "/", i, "_", comp, "_DIABLO_var_", keepvar, "_min.txt", sep=""
      )
      print("Writing DIABLO loadings to:")
      print(path_max)
      print(path_min)
      write.table(as.data.frame(loading_max),file=path_max,quote=FALSE,sep="\t")
      write.table(as.data.frame(loading_min),file=path_min,quote=FALSE,sep="\t")
    }
  }
}

assess_performance = function(data, dist, ncomp) {
  # review performance of diablo
  # remember to use the same distance metric which had the max value!
  # print("Assessing performance...")
  # perf_diablo = perf(data, validation='loo', M=10, nrepeat=10, dist=dist)
  # perf.diablo  # lists the different outputs

  # Performance with Majority vote
  # print(perf_diablo$MajorityVote.error.rate)

  # ROC and AUC criteria are not particularly insightful in relation to the
  # performance evaluation of our methods, but can complement the analysis.
  print("Plotting ROC...")
  mapply(function(x) auroc(data, x), seq(ncomp))
  return(perf_diablo)
}

predict_diablo = function(data, test, classes) {
  # prepare test set data: here one block (proteins) is missing
  print("Predicting data on an external test set...")
  predict.diablo = predict(data, newdata = test)
  # the warning message will inform us that one block is missing
  #predict.diablo # list the different outputs
  print("Getting confusion matrix...")
  confusion_matrix = get.confusion_matrix(
    truth=classes, predicted=predict.diablo$WeightedVote$max.dist[,2])
  print(confusion_matrix)
  print(get.BER(confusion_matrix))
}
