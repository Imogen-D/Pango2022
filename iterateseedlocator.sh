#!/bin/sh

echo load modules
#module unload bioinfo/locator-4760945
#module unload system/Miniconda3

module load system/Miniconda3
module load bioinfo/locator-4760945


indv=$1 ; echo $1
prefix=$2 ; echo $2
res=$3 ; echo $3
loc=$4 ; echo $4
vcffile=$5
metadata=$6
indv_loc_dir=$7

#indv_loc_dir=$res/locator/${indv}.${prefix}00_lociindv
#vcffile=$indv_loc_dir/${indv}lociindv.vcf.gz
#metadata=$indv_loc_dir/metadata.txt
loc_dir=$res/locator/seeding.concat.allWCA



#normal scirpt
while read seed ; do
  echo --vcf $vcffile --sample_data $metadata --out $indv_loc_dir/${indv}.${prefix}00_${seed}_lociindv
  x=$(ls $indv_loc_dir/*predlocs* | grep $seed)
  if [ -z "$x" ]
    then
    echo $seed
    locator.py --vcf $vcffile --sample_data $metadata --out $indv_loc_dir/${indv}.${prefix}00_${seed}_locindv --seed $seed
    rm $indv_loc_dir/*weights.hdf5
    else
    echo $seed no locator required
  fi
done < $loc_dir/seednumbers.txt

exit


#this loop is to fix those that were cancelled halfway, so see how far done and then iterate to suit
num=$(ls $indv_loc_dir/${indv}.${prefix}00_*_locindv*predlocs.txt | wc -l)
n=$((num+1))

tail -n "+$n" $loc_dir/seednumbers.txt > $indv_loc_dir/${indv}.${prefix}.tempseednumbers.txt

while IFS= read -r seed; do
  echo $seed
  echo --vcf $vcffile --sample_data $metadata --out $indv_loc_dir/${indv}.${prefix}00_${seed}_lociindv
  locator.py --vcf $vcffile --sample_data $metadata --out $indv_loc_dir/${indv}.${prefix}00_${seed}_locindv --seed $seed
done < $indv_loc_dir/${indv}.${prefix}.tempseednumbers.txt

rm $indv_loc_dir/tempseednumbers.txt

exit
