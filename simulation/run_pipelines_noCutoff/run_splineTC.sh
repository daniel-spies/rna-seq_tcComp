#!/bin/bash
dir=~/Documents/phd/data/simulation/study
dataDir=${dir}/run_data
scriptDir=~/Documents/phd/programming/scripts/RNAseq/simulation/study/run_pipelines_noCutoff

files=$(find $dataDir -type f -name "*.txt" | grep "contr" | grep -v "old" | grep -v "GSE")

mkdir -p ${dir}/results_noCutoff/splineTC
Rscript --vanilla $scriptDir/splineTC.R $files ${dir}/results_noCutoff/splineTC