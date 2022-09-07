#!/bin/sh

vcffile=$1 ; echo $1
loc_dir=$2 ; echo $2
loc=$3 ; echo $3
prefix=$4 ; echo $4
basedir=$5 ; echo $5
mm=$6 ; echo $6
samp=$7 ; echo $7

#loading modules
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/bcftools-1.3.1

#echo filtering 75 missing min mac 3

if  [[ ! -z "$mm" ]] ; then
  echo filtering all below ${mm}%
  vcftools --gzvcf ${vcffile} --remove-indels --min-alleles 2 --max-alleles 2  --max-missing 0.${mm} --mac 3 --recode --recode-INFO-all --stdout | bgzip -f -c > $loc_dir/${loc}${prefix}${mm}filtered1.vcf.gz
  echo tabix index with bcftools
  bcftools index -t -f $loc_dir/${loc}${prefix}${mm}filtered1.vcf.gz 
  echo cutting samples
  bcftools view  --force-samples -S $basedir/${prefix}.txt  $loc_dir/${loc}${prefix}${mm}filtered1.vcf.gz  | bgzip -f -c > ${loc_dir}/${loc}${prefix}${mm}filtered2.vcf.gz
  bcftools stats ${loc_dir}/${loc}${prefix}${mm}filtered2.vcf.gz > $loc_dir/${loc}${prefix}${mm}vcffiltered.stats
  mv ${loc_dir}/${loc}${prefix}${mm}filtered2.vcf.gz ${loc_dir}/${loc}${prefix}${mm}filtered.vcf.gz
  rm $loc_dir/*filtered1.vcf.gz
elif  [ -z "$mm" ] ; then
  echo not filtering
  vcftools --gzvcf ${vcffile}  --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | bgzip -f -c > $loc_dir/${loc}${prefix}${mm}filtered1.vcf.gz
  #bgzip -f -c $vcffile > $loc_dir/${loc}${prefix}filtered1.vcf.gz
  echo tabix index with bcftools
  bcftools index -t -f $loc_dir/${loc}${prefix}${mm}filtered1.vcf.gz #${vcffile} #
  echo cutting samples
  bcftools view --force-samples -S $basedir/${prefix}.txt  $loc_dir/${loc}${prefix}${mm}filtered1.vcf.gz   | bgzip -f -c > ${loc_dir}/${loc}${prefix}${mm}filtered2.vcf.gz 
  bcftools stats ${loc_dir}/${loc}${prefix}${mm}filtered2.vcf.gz > $loc_dir/${loc}${prefix}${samp}${mm}vcffiltered.stats
  mv ${loc_dir}/${loc}${prefix}${mm}filtered2.vcf.gz ${loc_dir}/${loc}${prefix}${samp}${mm}filtered.vcf.gz
  rm $loc_dir/*filtered1.vcf.gz
fi

exit
#WAS 0.75 only changed to 0.90 FOR MM90

echo tabix index with bcftools
bcftools index -t -f $loc_dir/${loc}${prefix}filtered1.vcf.gz

echo filtering on sample list $prefix and creating stat file
if [ $prefix = "all" ]; then #stops removal of any samples
  bcftools stats ${loc_dir}/${loc}${prefix}filtered1.vcf.gz > $loc_dir/${loc}${prefix}vcffiltered.stats
  mv ${loc_dir}/${loc}${prefix}filtered1.vcf.gz ${loc_dir}/${loc}${prefix}filtered.vcf.gz
elif [ $prefix != "all" ]; then
  bcftools view -S $basedir/${prefix}.txt  $loc_dir/${loc}${prefix}filtered1.vcf.gz  | bgzip -f -c > ${loc_dir}/${loc}${prefix}filtered2.vcf.gz
  bcftools stats ${loc_dir}/${loc}${prefix}filtered2.vcf.gz > $loc_dir/${loc}${prefix}vcffiltered.stats
  mv ${loc_dir}/${loc}${prefix}filtered2.vcf.gz ${loc_dir}/${loc}${prefix}filtered.vcf.gz
fi


exit

echo filtering for quality
#vcftools --gzvcf ${vcffile}.gz --remove-indels --minQ 30 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | gzip -c > ${vcffile}filtered.gz
#maf 0.002 is every alle 2 times; minimum depth of 3 reads on genotype 3 --max-missing 0.9  too stringent

echo filtering for mac
#vcftools --gzvcf ${vcffile}.gz --remove-indels --min-alleles 2 --max-alleles 2 --mac 3 --recode --recode-INFO-all --stdout | gzip -c > ${vcffile}filtered.gz

echo filtering for biallelic
#vcftools --gzvcf ${vcffile}.gz --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | gzip -c > ${vcffile}filtered.gz
#this doens't filter any out

echo filtering for all
#vcftools --gzvcf ${vcffile}.gz --remove-indels --maf 0.002 --minQ 30 --min-alleles 2 --max-alleles 2 --mac 3 --minDP 3 --max-missing 0.9 --recode --recode-INFO-all --stdout | gzip -c > ${vcffile}filtered.gz
#nothing left

echo filtering for all but quality and max missing # --max-missing 0.9
#vcftools --gzvcf ${vcffile}.gz --remove-indels --maf 0.002 --min-alleles 2 --max-alleles 2 --mac 3 --minDP 3 --recode --recode-INFO-all --stdout | gzip -c > ${vcffile}filtered.gz

echo tabix index with bcftools
bcftools index -t -f ${vcffile}filtered.gz

bcftools stats ${vcffile}filtered.gz > $loc_dir/vcffiltered.stats
exit


echo number lines before sorting
wc -l  $vcffile

echo sorting file #might need todo this sort -k1,1V -k2,2n my.vcf
#cat $vcffile | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > $vcffile
vcf-sort $vcffile > $vcffile.sorted
echo number lines once sorted
wc -l  $vcffile 
wc -l  $vcffile.sorted
mv  $vcffile.sorted $vcffile

echo zipping
bgzip -c $vcffile > $vcffile.gz

echo indexing
echo bcftools index
bcftools index -f ${vcffile}.gz
echo tabix index with bcftools
bcftools index -t -f ${vcffile}.gz

echo stat file on full file to $loc_dir/vcf.stats
bcftools stats $vcffile.gz > $loc_dir/vcf.stats



#this function isn't available so not run
echo missing indv
vcftools --gzvcf filtered_1_${vcffile} --missing_indv

##looking at indivudlas missing data from https://www.ddocent.com/filtering/
cat out.imiss

mawk '!/IN/' out.imiss | cut -f5 > totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
min reads 155000 min cov 1.5

mawk '$5 > 0.5' out.imiss | cut -f1 > ${loc}lowDP.indv

#can refilter with --min-meanDP and higher --max-missing with wanted
