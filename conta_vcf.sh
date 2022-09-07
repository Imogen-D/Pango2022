#!/bin/sh
#SBATCH -p workq


ddir=$1
id=$2

module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15

echo "run vcftools"
  echo $id
  time vcftools --vcf $ddir/populations.snps.vcf --recode --out $ddir/pop.$id --max-missing 0.95 --indv $id --minDP 2 --non-ref-af 0.05
echo "zip final vcf file"  
  time gzip $ddir/pop.$id.recode.vcf
echo "done"



