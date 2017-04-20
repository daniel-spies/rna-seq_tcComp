#!/bin/bash
# submit dummy job for queue name
bsub -J sim_study "sleep 3;echo "waited 3 sec""
# submit DyNB jobs
for file in $(ls ~/data)
do
    bsub -J sim_study -n1 -N -W 168:00 -w "ended(sim_study)" matlab -singleCompThread -nodisplay -r "DyNB_caller $file; quit"
done
# submit nsgp jobs
for file in $(ls ~/data)
do
    bsub -J sim_study -n1 -N -W 168:00 -w "ended(sim_study)" matlab -singleCompThread -nodisplay -r "nsgp_hmc $file;quit"
done

for file in $(ls ~/data)
do
    bsub -J sim_study -n1 -N -W 168:00 -w "ended(sim_study)" matlab -singleCompThread -nodisplay -r "nsgp_hmc $file;quit"
done

## test runs
folder=/cluster/home04/biol/dspies/data/simulation_study/30M_3rep_8TP
for idx in $(seq 1 10)
do
    file=30M_3rep_8TP_block_${idx}_formated.mat
    bsub -J sim_nsgp_hmc -n1 -N -W 168:00 matlab -singleCompThread -nodisplay -r "nsgp_hmc $folder/$file 16;quit"
done