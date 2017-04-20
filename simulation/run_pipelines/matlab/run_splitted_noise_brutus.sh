#!/bin/bash
# submitting files to brutus that have been splitted into subfiles
# each splitted fil list is run in nBins parallel queues
# next file will only start after the previous one has finished

# set splitting and running parameters
# the duration of each submitted job is calucated depending on the number of genes of each file
# assuming that 4 genes will be processed per hour
nCPU=16
walltime=16:00
nSplit=100
nBins=10
declare -a members=($(seq 0 $nBins $(($nSplit+1))))

previous_job=initial
bsub -J "initial" "sleep 3;echo "waited 3 sec""
for tool in DyNB_caller nsgp_hmc nsgp_lso
do
    for exp in $(find /cluster/home04/biol/dspies/data/simulation_study/noise -mindepth 1 -type d)
    do
        queue=$(basename $exp)
        for line in $(seq 1 $(expr $nSplit / $nBins))
        do
            for fileId in $(seq $((${members[$(($line-1))]}+1)) ${members[$line]})
            do
                file=$exp/${queue}_block_${fileId}_formated.mat
                bsub -J ${queue}_${line} -w "ended($previous_job)" -n 1 -N -W 120:00 matlab -singleCompThread -nodisplay -r "$tool $file $nCPU $walltime; quit"
                sleep 0.1 ## in order that jobs are registered as jobs with a unique id
            done
            previous_job=${queue}_${line}
            sleep 2 ## in order that jobs are already submited, otherwise the dependency might throw an error 
        done
    done
done