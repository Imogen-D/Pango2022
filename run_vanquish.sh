#!/bin/sh
#SBATCH -p workq

vdir=$1
rads=$2
DATA=$3
sample_file=$4
BIN=$5
x=$6

echo "run vanquish"
module load system/R-3.4.3 ; Rscript $BIN/vanquish_ploidy_plot.R $vdir $rads $sample_file $x

echo "done"

