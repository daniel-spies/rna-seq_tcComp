#!/bin/bash
cd ~/Documents/phd/data/simulation/study/results/timeSeq/

# merge result files
for folder in $(find * -type d)
do
    for file in $(ls $folder/*_timeSeq.txt)
    do
        cut -d" " -f2 $file | tail -n +2
    done > ${folder}_timeSeq.txt
done