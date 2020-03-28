#!/usr/bin/Rscript
# combine translatome and proteomics data for sars-cov-2
# data originally from DOI:10.21203/rs.3.rs-17218/v1 - supp tables 1 and 2
library(argparser, quietly=TRUE)
library(mixOmics)

parse_data = function(infile_path, offset=0) {
  # load in omics data into a diablo-compatible format
  print("Parsing file:")
  print(infile_path)
  return(read.table(infile_path, sep="\t", header=TRUE, row.names=1) + offset)
}

remove_novar = function(data) {
  # samples with zero variance are meaningless for PCA
  return(data[, which(apply(data, 2, var) != 0)])
}

parse_classes = function(infile_path) {
  # load in class data for diablo
  data = read.table(infile_path, sep="\t", header=TRUE, row.names=1)
  return(unlist(as.vector(t(data))))
}

create_design = function(data) {
  design = matrix(0.1, ncol = length(data), nrow = length(data),
                  dimnames = list(names(data), names(data)))
  diag(design) = 0
  return(design)
}

plot_individual_blocks = function(data, classes) {
  # do pca on individual classes before proceeding
  names = names(data)

  print("Removing 0 variance columns from data...")
  data = lapply(data, remove_novar)

  print("Showing PCA component contribution...")
  data_pca = lapply(data, pca, ncomp=dim(classes)[1], center=TRUE, scale=TRUE)

  print("Plotting PCA component contribution...")
  # lapply(data_pca, plot, title="Screeplot")
  mapply(function(x, y) plot(x, main=paste(y, "Screeplot")), data_pca, names)

  print("Plotting PCA...")
  # lapply(data_pca, plotIndiv, comp=c(1,2), ind.names=TRUE, group=classes,
  #   legend=TRUE, title="PCA 1/2")
  mapply(function(x, y) plotIndiv(x, comp=c(1,2), ind.names=TRUE, group=classes,
    legend=TRUE, title=paste(y, "PCA 1/2")), data_pca, names)

  print("Plotting correlation circle plots...")
  # lapply(data_pca, plotVar, comp=c(1, 2), var.names=TRUE, title="PCA 1/2")
  mapply(function(x, y) plotVar(x, comp=c(1, 2), title=paste(y, "PCA 1/2")),
    data_pca, names)

  print("Plotting biplots...")
  mapply(function(x, y, z) biplot(y, cex=0.7, xlabs=paste(classes, 1:nrow(x)),
    main=paste(z, "Biplot")), data, data_pca, names)
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
  test_keepX = list(proteome = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                    translatome = c(5:9, seq(10, 18, 2), seq(20,30,5)))
  tune_data = tune.block.splsda(X = data, Y = classes, ncomp = ncomp,
                                test.keepX = test_keepX, design = design,
                                validation = 'Mfold', folds = 10, nrepeat = 1,
                                cpus = cpus, dist = dist)
  list_keepX = tune_data$choice.keepX
  return(list_keepX)
}

run_diablo = function(data, classes, ncomp, keepx, design) {
  print("Running DIABLO...")
  sgccda_res = block.splsda(X = data, Y = classes, ncomp = ncomp,
                            keepX = keepx, design = design)
  return(sgccda_res)
}

plot_diablo = function(data) {
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
  # remember to use the same distance metric which had the max value
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
