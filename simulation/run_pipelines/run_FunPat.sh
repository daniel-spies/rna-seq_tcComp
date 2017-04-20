#!/bin/bash
dir=~/Documents/phd/data/simulation/study
cd $dir
dataDir=${dir}/run_data
scriptDir=~/Documents/phd/programming/scripts/RNAseq/simulation/study/run_pipelines
outDir=${dir}/results/FunPat

mkdir -p $outDir
files=$(ls $dataDir/*contr.txt | grep -v GSE)
Rscript --vanilla $scriptDir/FunPat.R $files $outDir

# clear data
rm -rf $dir/Gene_Ranking