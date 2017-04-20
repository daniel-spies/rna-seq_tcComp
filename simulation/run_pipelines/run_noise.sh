#!/bin/bash
path=~/Documents/phd/data/simulation/study
tools=$scripts/RNAseq/simulation/study/run_pipelines
outDir=$path/results/noise
files=$(ls $path/countTables/noise/* | grep contr)

## run TC tools on noisy data
mkdir -p $outDir/splineTC
Rscript $tools/splineTC.R $files $outDir/splineTC

mkdir -p $outDir/maSigPro
Rscript $tools/maSigPro.R $files $outDir/maSigPro

mkdir -p $outDir/EBSeqHMM
Rscript $tools/EBSeqHMM.R $files $outDir/EBSeqHMM

mkdir -p $outDir/pairwise
Rscript $tools/pairwise.R $files $outDir/pairwise

mkdir -p $outDir/FunPat
Rscript $tools/FunPat.R $files $outDir/FunPat