#!/bin/bash
dir=~/Documents/phd/data/simulation/study
dataDir=${dir}/run_data
scriptDir=~/Documents/phd/programming/scripts/RNAseq/simulation/study/run_pipelines

# standard is present in all simulation categories, filter out and run only once
files=$(find $dataDir -type f -name "*.txt" | grep "contr" )
Rscript --vanilla $scriptDir/pairwise.R $files ${dir}/results/pairwise