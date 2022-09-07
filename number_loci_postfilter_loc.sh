#!/bin/sh

samples=$1 ; echo $1
base_loc_dir=$2 ; echo $2
prefix=$3 ; echo $3
loc=$4 ; echo $4
res=$5 ; echo $5

rm $base_loc_dir/number_loci_individual_filtering.txt

while read indv ; do
VARS=(75 80 85 90)
for mnum in "${VARS[@]}"; do
  mm=mm${mnum}
  loc_dir=$res/locator/${loc}.${prefix}${mm}
  vcffile=${loc_dir}/${loc}${prefix}${mnum}filtered.vcf.gz
num=$(gunzip -c $vcffile | head -n 50 | grep CHROM |  awk -v b="$indv" '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}') #getting column number of infivual
#getting list of positions only found in that indivual (removing ./.)
x=$(gunzip -c $vcffile | awk -v c1=$num '{print $1, $2, $c1}' | sed 's/ /\t/g' | grep -v -F "./." | awk '{print $1, $2}' | sed 's/ /\t/g' | sed '1,29d' | wc -l) #numbe rof loci
printf "$indv\t$x\t$mnum\n" >> $base_loc_dir/number_loci_individual_filtering.txt
done
mnum=00
vcffile=$base_loc_dir/concat.allWCA/concatallWCAfiltered.vcf.gz
num=$(gunzip -c $vcffile | head -n 50 | grep CHROM |  awk -v b="$indv" '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}')
x=$(gunzip -c $vcffile | awk -v c1=$num '{print $1, $2, $c1}' | sed 's/ /\t/g' | grep -v -F "./." | awk '{print $1, $2}' | sed 's/ /\t/g' | sed '1,29d' | wc -l) #numbe rof loci
printf "$indv\t$x\t$mnum\n" >> $base_loc_dir/number_loci_individual_filtering.txt
done < $samples