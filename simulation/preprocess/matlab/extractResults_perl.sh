#!/bin/bash

cd ~/Documents/phd/data/simulation/study/results/noise/DyNB
for file in $(ls *.txt)
do
    awk 'function abs(x){return ((x < 0.0) ? -x : x)}{print abs(log($2)),$1}' $file | perl -e 'print sort { $a<=>$b } <>' | tail -n 1200 | awk '{print $2,$1}' > ${file%.txt}_result.txt
done

ls *.txt | grep -v result | xargs rm