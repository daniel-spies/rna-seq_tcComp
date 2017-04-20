#!/bin/bash
cd ~/Documents/phd/data/simulation/study/stats

# percentage of unique TP per tool
cat methods_overlap_TP.txt | tail -n +2 | cut -d" "  -f2-9 | 
    awk '{for(i=1;i<=NF;i++) t+=$i; print $NF/t; t=0}' | awk '{s+=$1} END {print s/NR}'
cat methods_overlap_TP_top5.txt | tail -n +2 | cut -d" "  -f2-9 | 
    awk '{for(i=1;i<=NF;i++) t+=$i; print $NF/t; t=0}' | awk '{s+=$1} END {print s/NR}'

# percentage of unique FP per tool
cat methods_overlap_FP.txt | tail -n +2 |head -n 5 | cut -d" "  -f2-9 | 
    awk '{for(i=1;i<=NF;i++) t+=$i; if(t > 0) print $NF/t; t=0}' | awk '{s+=$1} END {print s/NR}'
cat methods_overlap_FP_top5.txt | tail -n +2 | cut -d" "  -f2-9 | 
    awk '{for(i=1;i<=NF;i++) t+=$i; if(t > 0) print $NF/t; t=0}' | awk '{s+=$1} END {print s/NR}'