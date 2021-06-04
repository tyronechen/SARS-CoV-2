Rscript  inst/scripts/run_pipeline.R --tune_off \
  --classes ../data/case_study_1/classes_diablo.txt \
  --classes_secondary ../data/case_study_1/pch.txt \
  --dropna_classes TRUE \
  --dropna_prop 0.3 \
  --data ../data/case_study_1/diablo_proteome.txt ../data/case_study_1/diablo_translatome.txt \
  --data_names proteome translatome \
  --force_unique TRUE \
  --mappings ../data/case_study_1/proteome_mapfile.txt ../data/case_study_1/translatome_mapfile.txt \
  --ncpus 2 \
  --diablocomp 3 \
  --linkage 0.1 \
  --diablo_keepx 10 25 50 \
  --icomp 0 \
  --zero_as_na TRUE \
  --replace_missing FALSE \
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
  --outfile_dir ../results/case_study_1 \
  --corr_cutoff 0.6 \
  --rdata casestudy1.RData \
  --plot casestudy1.pdf \
  --args Rscript_cs1.sh
