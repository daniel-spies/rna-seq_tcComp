path=~/Documents/phd/data/simulation/study/results
script=~/Documents/phd/programming/scripts/RNAseq/simulation/study/
wd=$(pwd)

## extracting the top 2500 genes (abs(log) bayes factor) for matlab implementations
## opted for top 2500 as on average the other tools identfied 2460 candidates
cd $script
for method in dynb nsgp
do
    for file in $(ls $path/${method}/*.mat)
    do
        matlab -nosplash -nodesktop -nodisplay -nojvm -r "extractResults $file,quit"
        cat ${file%.mat}.txt | grep -v NaN | awk 'function abs(x){return ((x < 0.0) ? -x : x)}{print abs(log($2)),$1}' | 
            perl -e 'print sort { $a<=>$b } <>' | tail -n 1200 | awk '{print $2,$1}' > ${file%_result.mat}.txt
        rm ${file%.mat}.txt
    done
done
cd $wd