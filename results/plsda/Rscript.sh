#!/bin/bash
cd ../src
Rscript diablo_train.R \
   ../data/classes_diablo.txt  \
   --classes_secondary ../data/pch.txt \
   --dropna_classes TRUE \
   --dropna_prop 0 \
   --data ../data/diablo_proteome.txt ../data/diablo_translatome.txt \
   --mappings ../data/proteome_mapfile.txt ../data/translatome_mapfile.txt \
   --ncpus 6 \
   --dcomp 0 \
   --icomp 24 \
   --pcomp 10 \
   --plsdacomp 4 \
   --splsdacomp 4 \
   --splsda_keepx 10 50 100 250 \
   --mdist centroids \
   --contrib max \
   --outfile_dir ../results/centroids_i24_comp4_loadingmax \
   --rdata RData.RData \
   --plot Rplots.pdf \
   --args Rscript.sh
