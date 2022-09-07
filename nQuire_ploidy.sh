#!/bin/sh
#SBATCH -p workq

echo "################################################"
echo "assess individual ploidy data"
echo "################################################"

echo "################################################"
echo "1_loading modules"
module load bioinfo/nQuire-a990a88

echo "################################################"
echo "2 loading parameters"
nt=$1
data_dir=$2
indv=$3
sb_dir=$4
ref=$5
nquire_dir=$6
prefix=$7
suffix=$8
reads=$9
memo=${10}
memi=$(( memo / nt ))
#noMTnoCP=${11}

out_ext=.sorted_filtered_autosomes.bam ; #define the extensions of our files
#AN_123_2018_rmdups_noMTnoCP_spini_sorted.bam


echo "################################################"
echo "create .bin file from bam file"
nQuire create -b $sb_dir/${indv}${out_ext} -o $nquire_dir/${indv}_${prefix}.nquire -c 10 -q 60

echo "################################################"
echo "create histogram from .bin file"
nQuire histo $nquire_dir/${indv}_${prefix}.nquire.bin

echo "################################################"
echo "denoise .bin file"
#nQuire denoise $nquire_dir/${indv}_${prefix}.nquire.bin -o $nquire_dir/${indv}_${prefix}.nquire.denoised

echo "################################################"
echo "Assessing ploidy level with lrd model"
nQuire lrdmodel $nquire_dir/${indv}_${prefix}.nquire.bin

echo "################################################"
echo "Assessing ploidy level with modeltest"
nQuire modeltest $nquire_dir/${indv}_${prefix}.nquire.bin

echo "################################################"
echo "Assessing ploidy level with estmodel"
nQuire estmodel $nquire_dir/${indv}_${prefix}.nquire.bin

echo "################################################"
echo "Assessing ploidy level with histotest"
nQuire histotest $nquire_dir/${indv}_${prefix}.nquire.bin