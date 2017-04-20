#!/bin/bash
dir=~/Documents/phd/data/simulation/study
dataDir=${dir}/countTables
scriptDir=~/Documents/phd/programming/scripts/RNAseq/simulation/study/run_pipelines
standard=30M_3rep_4TP

files=$(find $dataDir -type f -name "*.txt" ! -path "$dataDir/standard/*" ! -path "$dataDir/raw/*" | grep "contr" | grep -v $standard)
file=$(echo $files | tr ' ' '\n' | grep -v 10M)
# user bigger stack and cell sizes to handle 12 TP data sets
Rscript --vanilla --max-ppsize=500000 --min-nsize=500000 --min-vsize=900000 $scriptDir/EBSeqHMM.R $file ${dir}/results/EBSeqHMM
Rscript --vanilla $scriptDir/EBSeqHMM.R $dataDir/standard/${standard}_contr.txt ${dir}/results/EBSeqHMM