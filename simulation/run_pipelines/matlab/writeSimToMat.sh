#!/bin/bash
dir=~/Documents/phd/data/simulation/study/countTables
out=~/Documents/phd/data/simulation/study/matlabTables
script=~/Documents/phd/programming/scripts/RNAseq/simulation/study/run_pipelines/matlab

Rscript --vanilla $script/writeSimToMat.R $dir $out