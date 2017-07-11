!/bin/bash
dir=/cclab_nas/dspies/data/simulated/study
dataDir=${dir}/run_data
scriptDir=${dir}/scripts

# standard scenarios
${dir}/results/lmms
files=$(find $dataDir -type f -name "*contr.txt" |  grep -v GSE | grep -v noise)
Rscript --vanilla $scriptDir/lmms.R $files ${dir}/results/lmms

# noise
mkdir ${dir}/results/noise/lmms
files=$(find $dataDir/noise -type f -name "*contr.txt")
Rscript --vanilla $scriptDir/lmms.R $files ${dir}/results/noise/lmms