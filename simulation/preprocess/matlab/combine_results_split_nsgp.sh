script=~/Documents/phd/programming/scripts/RNAseq/simulation/study
folder=~/Documents/phd/data/simulation/study/results/nsgp/30M_3rep_12TP
outFile=${folder}_nsgp_result
nBlocks=20
wd=$(pwd)

cd $script
matlab -nodisplay -r "combine_results_split_nsgp $folder $outFile $nBlocks;quit"

cd $wd