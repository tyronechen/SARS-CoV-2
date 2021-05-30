Rscript  inst/scripts/run_pipeline.R --tune_off \
  --classes ../data/case_study_2/classes_diablo.tsv \
  --dropna_classes FALSE \
  --dropna_prop 0.3 \
  --data \
    ../data/case_study_2/data_lipidomics.tsv \
    ../data/case_study_2/data_metabolomics.tsv \
    ../data/case_study_2/data_proteomics.tsv \
    ../data/case_study_2/data_transcriptomics_imputed_all.tsv \
  --data_names lipidome metabolome proteome transcriptome \
  --mappings NA \
  --ncpus 2 \
  --diablocomp 3 \
  --linkage 0.5 \
  --diablo_keepx 10 25 50 \
  --icomp 0 \
  --pcomp 3 \
  --plsdacomp 3 \
  --splsdacomp 3 \
  --splsda_keepx 10 25 50 \
  --dist_plsda centroids.dist \
  --dist_splsda centroids.dist \
  --dist_diablo mahalanobis.dist \
  --cross_val loo \
  --cross_val_nrepeat 1 \
  --cross_val_folds 10 \
  --contrib max \
  --outfile_dir ../results/case_study_2 \
  --corr_cutoff 0.6 \
  --rdata casestudy2.RData \
  --plot casestudy2.pdf \
  --args Rscript_cs2.sh
