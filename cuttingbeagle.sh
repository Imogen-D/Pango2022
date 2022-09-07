#!/bin/sh
#SBATCH -p workq


loc=$1
res=$2
ngsadmix_dir=$3
output=$4
input=$5

echo "remove low coverage and lineages from beagle file"
module load system/R-3.4.3 ; Rscript /work/idumville/pangolins/RAD/bin/cuttingbeagle.R ${loc} $output $input $ngsadmix_dir/unfilteredall${loc}.beagle.gz $res/ngsadmix.autosomes/sampledataautosomes.txt $ngsadmix_dir
echo "done"