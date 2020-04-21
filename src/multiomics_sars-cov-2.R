#!/usr/bin/Rscript
# combine translatome and proteomics data for sars-cov-2
# data originally from DOI:10.21203/rs.3.rs-17218/v1 - supp tables 1 and 2
library(argparser, quietly=TRUE)
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

show_na_prop = function(data_na, name) {
  # first, sum the number of missing values per variable
  sum_na_per_var = apply(data_na, 2, function(x) {sum(is.na(x))})
  ​
  # simple plot shows that some variable have an NA rate greater than 30%
  plot(sum_na_per_var/ncol(data_na), type='h', xlab='variable index',
    ylab='NA rate', main=paste(name, 'NA rate per variable on unfiltered data'))
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

create_design = function(data) {
  # create design matrix from data: dataframe -> matrix (design)
  design = matrix(0.1, ncol = length(data), nrow = length(data),
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

impute_missing_ = function(data, ncomp=10, block_name="") {
  # impute missing values with nipals: dataframe (with NA) -> dataframe
  print("Number of components for imputation:")
  print(ncomp)
  nipals_tune = nipals(data, reconst=TRUE, ncomp=ncomp)
  barplot(nipals_tune$eig, main=paste(block_name, "Screeplot (nipals imputed)"),
    xlab="Number of components", ylab="Explained variance"
  )
  return(nipals_tune$rec)
}

impute_missing = function(data, ncomps) {
  mapply(function(x, y, z) impute_missing_(x, y, z), data, ncomps, names(data))
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
    mapply(function(x, y) plotIndiv(x, comp=c(1,2), ind.names=TRUE,
      group=classes, legend=TRUE, ncomp=ncomp,
      title=paste(title, y, "PCA 1/2"), pch=pch), data_pca, names)
    mapply(function(x, y) plotIndiv(x, comp=c(1,3), ind.names=TRUE,
      group=classes, legend=TRUE, ncomp=ncomp,
      title=paste(title, y, "PCA 1/3"), pch=pch), data_pca, names)
    mapply(function(x, y) plotIndiv(x, comp=c(2,3), ind.names=TRUE,
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
  mapply(function(x, y) plotVar(x, comp=c(1, 2), title=paste(y, "PCA 1/2")),
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
  mapply(function(x, y) plotIndiv(x, comp=c(1,2), ind.names=TRUE,
    group=classes, legend=TRUE, ncomp=ncomp,
    title=paste(title, y, "PCA M 1/2"), pch=pch), data_pca, names)
  mapply(function(x, y) plotIndiv(x, comp=c(1,3), ind.names=TRUE,
    group=classes, legend=TRUE, ncomp=ncomp,
    title=paste(title, y, "PCA M 1/3"), pch=pch), data_pca, names)
  mapply(function(x, y) plotIndiv(x, comp=c(2,3), ind.names=TRUE,
    group=classes, legend=TRUE, ncomp=ncomp,
    title=paste(title, y, "PCA M 2/3"), pch=pch), data_pca, names)
  return(data_pca)
}

plsda_classify = function(data, classes, pch=NA, title="", ncomp=0) {
  mapply(function(x, y) plsda_classify_(x, classes, pch, y, ncomp), data, title)
}

plsda_classify_ = function(data, classes, pch=NA, title="", ncomp=0) {
  # discriminate samples: list, vector, bool, integer -> list
  # single or multilevel PLS-DA
  if (!is.na(pch)) {
    print("Plotting multi level partial least squares discriminant analysis")
    data_plsda = plsda(data, Y=classes, multilevel=c(as.factor(pch)), ncomp=ncomp)
    plotIndiv(data_plsda, ind.names=TRUE, group=classes, legend=TRUE,
      pch=pch, title=paste(title, "PLSDA multi 1/2"), comp=c(1,2)
    )
    plotIndiv(data_plsda, ind.names=TRUE, group=classes, legend=TRUE,
      pch=pch, title=paste(title, "PLSDA multi 1/3"), comp=c(1,3)
    )
    plotIndiv(data_plsda, ind.names=TRUE, group=classes, legend=TRUE,
      pch=pch, title=paste(title, "PLSDA multi 2/3"), comp=c(2,3)
    )
  } else {
    print("Plotting single level partial least squares discriminant analysis")
    data_plsda = plsda(data, Y=classes, ncomp=ncomp)
    plotIndiv(data_plsda, ind.names=TRUE, group=classes, legend=TRUE,
      title=paste(title, "PLSDA single 1/2"), comp=c(1,2)
    )
    plotIndiv(data_plsda, ind.names=TRUE, group=classes, legend=TRUE,
      title=paste(title, "PLSDA single 1/3"), comp=c(1,3)
    )
    plotIndiv(data_plsda, ind.names=TRUE, group=classes, legend=TRUE,
      title=paste(title, "PLSDA single 2/3"), comp=c(2,3)
    )
  }
  mapply(function(x) auroc(data_plsda, roc.comp=x), seq(ncomp))
  loadings = plotLoadings(data_plsda,contrib='max',comp=ncomp,method='median',
    ndisplay=50, name.var=colnames(data), title=paste(title, "PLSDA loadings")
  )
  print(loadings)
  path = paste(title, "_PLSDA", ".txt", sep="")
  write.table(as.data.frame(loadings), file=path, quote=FALSE, sep="\t")
  return(data_plsda)
}

splsda_tune = function(data, classes, names, multilevel, ncomp=3, nrepeat=10,
  logratio="none", test_keepX=c(5,50,100), validation="loo", folds=10,
  dist="max.dist", cpus=2, progressBar=TRUE) {
    mapply(function(x, y) splsda_tune_(x, classes, names, multilevel, ncomp,
      nrepeat, logratio, test_keepX, validation, folds, dist, cpus, progressBar),
      data, names, SIMPLIFY=FALSE)
  }

splsda_tune_ = function(data, classes, names, multilevel, ncomp=3, nrepeat=10,
  logratio="none", test_keepX=c(5,50,100), validation="loo", folds=10,
  dist="max.dist", cpus=2, progressBar=TRUE) {
  # tune splsda components
  tuned = tune.splsda(data, Y=classes, multilevel=multilevel, ncomp=ncomp,
    nrepeat=nrepeat, logratio=logratio, test.keepX=test_keepX,
    validation=validation, folds=folds, dist=dist, cpus=cpus,
    progressBar=progressBar
  )
  print(plot(tuned, main=names))
  return(tuned)
}

splsda_classify = function(data, classes, pch=NA, title="", ncomp=0, keepX=NULL) {
  mapply(function(x, y, c, k) splsda_classify_(
    x, classes, pch, y, c, k
  ), data, title, ncomp, keepX)
}

splsda_classify_ = function(data, classes, pch=NA, title="", ncomp=NULL, keepX=NULL) {
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
  data_splsda = splsda(data, Y=classes, multilevel=pch, ncomp=ncomp, keepX=keepX)
  plotIndiv(data_splsda, ind.names=TRUE, group=classes, legend=TRUE,
    pch=pch, title=paste(title, "sPLSDA multi 1/2"), comp=c(1,2)
  )
  plotIndiv(data_splsda, ind.names=TRUE, group=classes, legend=TRUE,
    pch=pch, title=paste(title, "sPLSDA multi 1/3"), comp=c(1,3)
  )
  plotIndiv(data_splsda, ind.names=TRUE, group=classes, legend=TRUE,
    pch=pch, title=paste(title, "sPLSDA multi 2/3"), comp=c(2,3)
  )
  mapply(function(x) auroc(data_splsda, roc.comp=x), seq(ncomp))
  loadings = plotLoadings(data_splsda,contrib='max',comp=ncomp,method='median',
    ndisplay=50, name.var=colnames(data), title=paste(title, "sPLSDA loadings")
  )
  print(loadings)
  path = paste(title, "_sPLSDA", ".txt", sep="")
  write.table(as.data.frame(loadings), file=path, quote=FALSE, sep="\t")
  return(data_splsda)
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

tune_ncomp = function(data, classes, design) {
  # First, we fit a DIABLO model without variable selection to assess the global
  # performance and choose the number of components for the final DIABLO model.
  # The function perf is run with 10-fold cross validation repeated 10 times.
  # ncomp = length(unique(classes))
  print("Finding optimal number of components...")
  ncomp = 10
  sgccda_res = block.splsda(X = data, Y = classes, ncomp = ncomp, design = design)

  # this code takes a couple of min to run
  perf_diablo = perf(sgccda_res, validation = 'Mfold', folds = 10, nrepeat = 10)

  # print(perf.diablo)  # lists the different outputs
  plot(perf_diablo)
  # perf_diablo$choice.ncomp$WeightedVote
  print(perf_diablo$choice.ncomp)
  return(perf_diablo)
}

tune_keepx = function(data, classes, ncomp, design, cpus=2, dist="centroids.dist") {
  # This tuning function should be used to tune the keepX parameters in the
  #   block.splsda function.
  # We choose the optimal number of variables to select in each data set using
  # the tune function, for a grid of keepX values. Note that the function has
  # been set to favor the small-ish signature while allowing to obtain a
  # sufficient number of variables for downstream validation / interpretation.
  # See ?tune.block.splsda.
  print("Tuning keepX parameter...")
  # test_keepX = list(proteome = c(5:9, seq(10, 18, 2), seq(20,30,5)),
  #                   translatome = c(5:9, seq(10, 18, 2), seq(20,30,5)))
  test_keepX = mapply(function(name, dims) list(name=dims), names(data),
    rep(list(c(5:9, seq(10, 18, 2), seq(20,30,5))))
  )

  tune_data = tune.block.splsda(X = data, Y = classes, ncomp = ncomp,
                                test.keepX = test_keepX, design = design,
                                validation = 'Mfold', folds = 10, nrepeat = 1,
                                cpus = cpus, dist = dist)
  list_keepX = tune_data$choice.keepX
  return(list_keepX)
}

run_diablo = function(data, classes, ncomp, keepx, design) {
  # this is the actual part where diablo is run
  print("Running DIABLO...")
  sgccda_res = block.splsda(X = data, Y = classes, ncomp = ncomp,
                            keepX = keepx, design = design)
  return(sgccda_res)
}

plot_diablo = function(data) {
  # plot the diablo data with a series of diagnostic plots
  print("Plotting correlation betweem components...")
  plotDiablo(data, ncomp = 1)
  print("Plotting individual samples into space spanned by block components...")
  plotIndiv(data, ind.names = FALSE, legend = TRUE, title = 'DIABLO', ellipse = TRUE)
  print("Plotting arrow plot...")
  plotArrow(data, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
  print("Plotting correlation circle plot...")
  plotVar(data, var.names = FALSE, style = 'graphics', legend = TRUE,
    pch=c(16, 17), cex=c(2,2), col=c('darkorchid', 'lightgreen')
  )
  print("Plotting circos from similarity matrix...")
  circosPlot(data, cutoff = 0.7, line = TRUE,
             color.blocks= c('darkorchid', 'lightgreen'),
             color.cor = c("chocolate3","grey20"), size.labels = 1.5)
  print("Plotting relevance network from similarity matrix...")
  network(data, blocks = c(1,2),
          color.node = c('darkorchid', 'lightgreen'), cutoff = 0.4)
  print("Plotting loading weight of selected variables on each component and dataset...")
  plotLoadings(data, comp = 1, contrib = 'max', method = 'median')
  print("Plotting heatmap...")
  cimDiablo(data)
}

assess_performance = function(data, dist) {
  # review performance of diablo
  # remember to use the same distance metric which had the max value!
  print("Assessing performance...")
  perf_diablo = perf(data, validation='Mfold', M=10, nrepeat=10, dist=dist)
  #perf.diablo  # lists the different outputs

  # Performance with Majority vote
  # print(perf_diablo$MajorityVote.error.rate)

  # ROC and AUC criteria are not particularly insightful in relation to the
  # performance evaluation of our methods, but can complement the analysis.
  print("Plotting ROC...")
  auc_splsda = auroc(data, roc.block = "proteome", roc.comp = 1)
}

predict_diablo = function(data, test, classes) {
  # prepare test set data: here one block (proteins) is missing
  # data.test.TCGA = list(mRNA = breast.TCGA$data.test$mrna,
  #                       miRNA = breast.TCGA$data.test$mirna)
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
