#!/bin/sh
#SBATCH -p workq


beaglefile=$1
outfile=$2
nt=3
npc=$4
pcangsdir=$5
ib=$6
bin=$7
infofile=$8
loc=$loc

echo "load modules"
module load compiler/gcc-6.1.0
module load bioinfo/pcangsd-0.97

python /usr/local/bioinfo/src/PCAngsd/pcangsd-0.97/pcangsd.py -beagle $beaglefile -o $outfile -threads $nt -inbreed $ib

module load system/R-3.5.2 ; Rscript $bin/plot_pcangsd.R $infofile $pcangsdir $outfile $loc $npc

