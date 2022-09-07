#!/bin/sh

input=$1 ; echo input $1
num=$2 ; echo loci $2
suffix=$3 ; echo dataset $3
col=$4 #number of columns in dataset minus 4
filt=$5 #filtering
vcf=$6
minmac=$7

user=idumville
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin

as_dir=$res/assignment

module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/bcftools-1.3.1
module load system/datamash-1.3

#echo bgziping the original file
#gunzip -c $vcf.gz | bgzip -c > ${vcf}bgzip.gz

#echo indexing file
#tabix ${vcf}bgzip.gz
 

echo filtering
#vcftools --gzvcf ${vcf}bgzip.gz --max-missing 0.${filt} --recode --recode-INFO-all --stdout | bgzip -f -c > $as_dir/${filt}RUBIloci${minmac}.vcf.gz

echo indexing file
#tabix $as_dir/${filt}RUBIloci${minmac}.vcf.gz

echo getting samples
#bcftools query -l $as_dir/${filt}RUBIloci${minmac}.vcf.gz > $as_dir/rubiassamples${minmac}${filt}.txt

echo editing this file to add reference mixture etc
#rm $as_dir/rubiasmeta${minmac}${filt}.txt
#while read indv ; do
#  grep -w $indv $as_dir/rubiasmeta2concatfull.txt >> $as_dir/rubiasmeta${minmac}${filt}.txt
#done < $as_dir/rubiassamples${minmac}${filt}.txt

#rm $as_dir/rubiassamples${minmac}${filt}.txt

#sed -i '1 i\sample_type\trepunit\tcollection\tindiv' $as_dir/rubiasmeta${minmac}${filt}.txt

#echo Querying alleles
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' $as_dir/${filt}RUBIloci${minmac}.vcf.gz | sed 's/\t/_/1' | sed 's|\t|/c\t|1' > $as_dir/${filt}RUBIqueried${minmac}.vcf

echo changing alleles to nucleotides
#rm $as_dir/${filt}RUBIqueried${minmac}edit.vcf
#while read line; do
#  R=$(echo $line | awk '{print $2}')
#  A=$(echo $line | awk '{print $3}')
#  echo $line | sed -e "s|0/|${R}/|g" -e "s|1/|${A}/|g"  -e "s|/1|/${A}|g" -e "s|/0|/${R}|g" >> $as_dir/${filt}RUBIqueried${minmac}edit.vcf
#done < $as_dir/${filt}RUBIqueried${minmac}.vcf

echo remove intemediate
#rm $as_dir/${filt}RUBIqueried${minmac}.vcf

echo changing alleles to numbers, transposing
#sed 's|A/|1/|g;s|/A|/1|g;s|C/|2/|g;s|/C|/2|g;s|G/|3/|g;s|/G|/3|g;s|T/|4/|g;s|/T|/4|g' $as_dir/${filt}RUBIqueried${minmac}edit.vcf | cut -d' ' --complement -f2,3 | datamash transpose -t ' ' | sed 's/\./NA/g' | sed 's|/|\t|g' | sed 's/ /\t/g' > $as_dir/rubiasgeno${filt}${minmac}.txt

echo remove intemediate
#rm $as_dir/${filt}RUBIqueried${minmac}edit.vcf 

echo pasting the files
#paste $as_dir/rubiasmeta${minmac}${filt}.txt $as_dir/rubiasgeno${filt}${minmac}.txt | sed 's|/|\t|g' > $as_dir/rubiasinput${filt}${minmac}.txt

echo remove intemediates
#rm $as_dir/rubiasmeta${minmac}${filt}.txt 
#rm $as_dir/rubiasgeno${filt}${minmac}.txt

input=$as_dir/rubiasinput${filt}${minmac}.txt
col1=$(head -n 1 $input | awk --field-separator="\t" "{ print NF }")
col=$(expr $col1 - 4)

mkdir -p $as_dir/${filt}${minmac}assignments

echo now into R
module load system/R-4.1.2_gcc-9.3.0

Rscript /work/idumville/pangolins/RAD/bin/run.assignment.R $input ${filt} ${minmac} $col

exit 



