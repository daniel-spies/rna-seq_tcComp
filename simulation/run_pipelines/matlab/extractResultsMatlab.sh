path=~/Documents/phd/data/simulation/study/results_noCutoff
script=~/Documents/phd/programming/scripts/RNAseq/simulation/study/run_pipelines_noCutoff
wd=$(pwd)

cd $script
for method in dynb nsgp
do
    for file in $(ls $path/${method}/*.mat)
    do
        matlab -nosplash -nodesktop -nodisplay -nojvm -r "extractResultsMatlab $file,quit"
        cat ${file%.mat}.txt | grep -v NaN | awk 'function abs(x){return ((x < 0.0) ? -x : x)}{print abs(log($2)),$1}' | 
            perl -e 'print sort { $b<=>$a } <>' | awk '{print $2,$1}' > ${file%_result.mat}.txt
        rm ${file%.mat}.txt
    done
done
cd $wd