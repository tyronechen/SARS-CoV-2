Rscript diablo_train.R \
   ../data/classes_diablo.txt  \
   --args Rscript.sh \
   --classes_secondary ../data/pch.txt \
   --data ../data/diablo_proteome.txt ../data/diablo_translatome.txt \
   --dropna_classes TRUE \
   --dropna_prop 0 \
   --ncpus 6 \
   --dcomp 0 \
   --icomp 5 \
   --rdata RData.RData \
   --plot Rplots.pdf \
   --outfile_dir test \
   --pcomp 10 \
   --plsdacomp 4 \
   --splsdacomp 4 \
   --splsda_keepx 10,50,100,250 \
   --mdist mahalanobis \
   --contrib min
