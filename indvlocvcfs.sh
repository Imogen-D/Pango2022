#!/bin/sh

vcffile=$1 ; echo $1
indv_loc_dir=$2 ; echo $2
indv=$3 ; echo $3
oldmetadata=$4 ; echo $4
rurlist=$5 ; echo $5

#loading modules
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/bcftools-1.3.1

rm $indv_loc_dir/metadata.txt
touch $indv_loc_dir/metadata.txt

#making metadata
awk '{print $1}' $rurlist | while read rurindv; do 
x=$(grep -w $rurindv $oldmetadata)

printf "$x\n" | sed 's/ /\t/g'  >> $indv_loc_dir/metadata.txt
done
printf "$indv NA NA\n" | sed 's/ /\t/g'  >> $indv_loc_dir/metadata.txt
sed -i '1 i\sampleID\tx\ty'  $indv_loc_dir/metadata.txt
sed -i '/^$/d' $indv_loc_dir/metadata.txt

awk '{print $1}'  $indv_loc_dir/metadata.txt > $indv_loc_dir/samplesincl.txt

num=$(gunzip -c $vcffile | head -n 10000 | grep CHROM | awk -v b="$indv" '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}') #getting column number of infivual
#getting list of positions only found in that indivual (removing ./.)
gunzip -c $vcffile | awk -v c1=$num '{print $1, $2, $c1}' | sed 's/ /\t/g' | grep -v -F "./." | awk '{print $1, $2}' | sed 's/ /\t/g' > $indv_loc_dir/locilist.txt
#trimming vcf
vcftools --gzvcf ${vcffile} --positions $indv_loc_dir/locilist.txt --keep $indv_loc_dir/samplesincl.txt --recode --recode-INFO-all --stdout | bgzip -f -c > $indv_loc_dir/${indv}lociindv.vcf.gz

#rm $indv_loc_dir/samplesincl.txt