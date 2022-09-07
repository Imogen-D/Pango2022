#!/bin/sh

nt=$1 ; echo $1
popmap=$2 ; echo $2
stdir=$3 ; echo $3
popdir=$4 ; echo $4
r=$5
R=$6

echo loading modules
module load compiler/gcc-7.2.0
module load bioinfo/stacks-2.5

#populations -P $stdir/ -M $popmap -O $popdir -t $nt --min-mac 2 -p 1 -r $r -R $R --vcf --no-hap-exports

echo running populations
populations -P $stdir/ --popmap $popmap  -O $popdir -p 1  --min-mac 2 -r $r -R $R -t $nt --structure --genepop --write-single-snp  --no-hap-exports





#--plink 

#gzip $popdir/*vcf


#run again with STRUCTURE requirements // still need to delete first line once produced
#--structure
#--write_single_snp
#--ordered_export