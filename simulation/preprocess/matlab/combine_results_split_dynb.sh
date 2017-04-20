scriptDir=~/Documents/phd/programming/scripts/RNAseq/simulation/study
#dataDir=~/Documents/phd/data/simulation/study/results/noise/DyNB
dataDir=~/Documents/phd/data/simulation/study/results/dynb
nBlocks=20
wd=$(pwd)

cd $scriptDir
for folder in $(find $dataDir -type d -mindepth 1 -maxdepth 1 -exec basename {} \;)
do
    #outFile=$dataDir/${folder}_noise_dynb
    outFile=$dataDir/${folder}_dynb
    matlab -nodisplay -r "combine_results_split_dynb $dataDir/$folder $outFile $nBlocks;quit"
done
cd $wd