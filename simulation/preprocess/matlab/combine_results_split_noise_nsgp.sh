scriptDir=~/Documents/phd/programming/scripts/RNAseq/simulation/study
nBlocks=10
wd=$(pwd)

cd $scriptDir
for folder in $(find ~/Documents/phd/data/simulation/study/results/noise/nsgp -type d -mindepth 1)
do
    outFile=${folder}_noise_nsgp_result
    matlab -nodisplay -r "combine_results_split_nsgp $folder $outFile $nBlocks;quit"
done
cd $wd