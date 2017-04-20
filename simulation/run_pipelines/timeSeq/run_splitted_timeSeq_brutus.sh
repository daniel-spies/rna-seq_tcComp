#!/bin/bash
# run timeSeq splitted files
nCPU=16
walltime=24:00
previous_job=initial
bsub -J $previous_job "sleep 3;echo "waited 3 sec""

for exp in $(find timeSeq/data/* -type d)
do
    queue=$(basename $exp)
    for file in $(ls $exp)
    do
        bsub -J ${queue} -w "ended($previous_job)" -n $nCPU -N -W $walltime -R "rusage[mem=4096]" "Rscript --vanilla --slave timeSeq.R  timeSeq/data/$queue/$file"
    done
    previous_job=${queue}
done