#!/bin/bash
dir=~/Documents/phd/data/simulation/study/matlabTables
script=~~/Documents/phd/programming/scripts/RNAseq/simulation/study/matlab
matlab -nodisplay -r "try, $script/formatReadCounts('$dir'); catch, disp('failed'), end, quit"