#!/bin/bash

dir=~/Documents/phd/data/simulation/study/matlabTables
script=~/Documents/phd/programming/scripts/RNAseq/simulation/study/run_pipelines/matlab
split=100
wd=$(pwd)

cd $script
for folder in $(find $dir -type d -mindepth 1 -maxdepth 1| grep -v raw)
do
    matlab -nodisplay -r "formatReadCounts_split $folder $split;quit"
done
cd $wd