#!/bin/bash
# example usage of code
Rscript diablo_train.R ../data/classes_diablo.txt \
	--classes_secondary ../data/pch.txt \
	-d ../data/diablo_proteome.txt ../data/diablo_translatome.txt \
	--rdata RData.RData \
	-i 24 \
	-p Rplots.pdf \
	--dropna_classes TRUE \
	--pcomp 10 \
	--plsdacomp 10 \
	--splsdacomp 5 \
	--mdist mahalanobis \
	--ncpus 6 \
	--splsda_keepx 5 10 15 20 25
