#!/bin/bash
dir=/cclab_nas/dspies/data/simulated/study
dataDir=${dir}/run_data
scriptDir=${dir}/scripts
standard=30M_3rep_4TP

files=$(find $dataDir -type f -name "*contr.txt" | grep -v GSE)
# user bigger stack and cell sizes to handle 12 TP data sets
Rscript --vanilla --max-ppsize=500000 --min-nsize=500000 --min-vsize=900000 $scriptDir/EBSeqHMM.R $files ${dir}/results_noCutoff/EBSeqHMM
Rscript --vanilla $scriptDir/EBSeqHMM.R $dataDir/standard/${standard}_contr.txt ${dir}/results_noCutoff/EBSeqHMM