#!/bin/sh
Rscript -e 'library(multiomics); data(BPH2819); export <- function(name, data) {write.table(data.frame(data), paste(name, ".tsv", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)}; mapply(export, names(BPH2819), BPH2819, SIMPLIFY=FALSE)'

Rscript inst/scripts/run_pipeline.R \
  --classes classes.tsv \
  --data metabolome.tsv \
         proteome.tsv \
         transcriptome.tsv \
  --data_names metabolome proteome transcriptome \
  --ncpus 2 \
  --icomp 12 \
  --pcomp 10 \
  --plsdacomp 2 \
  --splsdacomp 2 \
  --diablocomp 2 \
  --dist_plsda "centroids.dist" \
  --dist_splsda "centroids.dist" \
  --dist_diablo "centroids.dist" \
  --cross_val "Mfold" \
  --cross_val_folds 5 \
  --cross_val_nrepeat 50 \
  --corr_cutoff 0.1 \
  --outfile_dir BPH2819 \
  --contrib "max" \
  --progress_bar
