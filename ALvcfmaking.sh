#!/bin/bash

vcffile=$1 ; echo $1
num=$2 ; echo $2
min=$3 ; echo $3
al_dir=$4 ; echo $4
bin=$5 ; echo $5

echo loading modules
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/bcftools-1.3.1

#vcftools --gzvcf ${vcffile} --keep $al_dir/allsamples.txt --max-missing 0.${num} --recode --recode-INFO-all --stdout | bgzip -f -c > $al_dir/${num}ALloci${min}.vcf.gz

#vcftools --gzvcf ${vcffile} --keep "/work/idumville/pangolins/RAD/data/WCAsamples.txt" --max-missing 0.${num} --recode --recode-INFO-all --stdout | bgzip -f -c > $al_dir/${num}WCAloci${min}.vcf.gz

echo locator
loc=concat
prefix=ALfilt${num}${min}
#prefix=WCAfilt${num}${min}
echo $prefix
VAR=all
$bin/iterateseedlocator.sh $VAR $prefix /work/idumville/pangolins/RAD/results $loc $al_dir/${num}ALloci${min}.vcf.gz $al_dir/metadata.txt $al_dir/filt${num}${min}/


#$bin/iterateseedlocator.sh $VAR $prefix /work/idumville/pangolins/RAD/results $loc $al_dir/${num}WCAloci${min}.vcf.gz $al_dir/WCAminmac2metadata.txt $al_dir
