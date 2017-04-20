#!/bin/bash

path=~/Documents/phd/data/simulation/study
tools=$scripts/RNAseq/simulation/study/run_pipelines
outDir=$path/results/GSE69822

## run TC tools on public data set with PTEN-KO / PIK3A mutation and within WT time series
for norm in $(ls $path/run_data/GSE69822/*_norm_contr.txt)
do
    file=${norm%_norm_contr.txt}_contr.txt
    exp=$(basename $norm | sed 's/.*_GSE69822_\(.*\)_norm_contr.txt/\1/')
    
    mkdir -p $outDir/$exp/splineTC
    Rscript $tools/splineTC.R $norm $outDir/$exp/splineTC

    #mkdir -p $outDir/$exp/maSigPro
    #Rscript $tools/maSigPro.R $norm $outDir/$exp/maSigPro

    #mkdir -p $outDir/$exp/FunPat
    #Rscript $tools/FunPat.R $file $outDir/$exp/FunPat

    #mkdir -p $outDir/$exp/EBSeqHMM
    #Rscript $tools/EBSeqHMM.R $file $outDir/$exp/EBSeqHMM

    mkdir -p $outDir/$exp/pairwise_DESeq2
    Rscript $tools/pairwise_DESeq2.R $file $outDir/$exp/pairwise_DESeq2

    mkdir -p $outDir/$exp/pairwise_edgeR
    Rscript $tools/pairwise.R $file $outDir/$exp/pairwise_edgeR
done