#!/bin/bash
dir=/cclab_nas/dspies/data/simulated/study
dataDir=${dir}/data
scriptDir=${dir}/scripts
standard=30M_3rep_4TP

# standard is present in all simulation categories, filter out and run only once
files=$(find $dataDir -type f -name "*.txt" ! -path "$dataDir/standard/*" ! -path "$dataDir/raw/*" | grep "contr" | grep -v $standard | grep -v GSE)
Rscript --vanilla $scriptDir/maSigPro_multiTimeCourse.R $files ${dir}/results/maSigPro &
Rscript --vanilla $scriptDir/maSigPro_multiTimeCourse.R $dataDir/standard/${standard}_contr.txt ${dir}/results/maSigPro &