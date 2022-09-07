#!/bin/sh


module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/bcftools-1.3.1

res=/work/idumville/pangolins/RAD/results

echo step 1
gunzip -c $res/locator/concat.allWCA/concatallWCAfiltered.vcf.gz | head -n 100 | grep "^#" | bgzip -c -f > $res/locator/concat.allWCA/concatallWCAfiltered2.vcf.gz
echo sorting to unzipped
gunzip -c $res/locator/concat.allWCA/concatallWCAfiltered.vcf.gz | grep -v "^#" | sort -k1,1V -k2,2g >> $res/locator/concat.allWCA/concatallWCAfiltered2temp.vcf

echo bgzipping and appending
bgzip -c -f $res/locator/concat.allWCA/concatallWCAfiltered2temp.vcf >> $res/locator/concat.allWCA/concatallWCAfiltered2.vcf.gz

echo removing unzipped 
rm $res/locator/concat.allWCA/concatallWCAfiltered2temp.vcf

echo indexing
bcftools index -t -f $res/locator/concat.allWCA/concatallWCAfiltered2.vcf.gz


tabix -f -p vcf $res/locator/concat.allWCA/concatallWCAfiltered2.vcf.gz

tabix -p vcf $res/locator/concatWCAUlowread.vcf.gz

echo merging #neither of these are working
bcftools merge $res/locator/concat.allWCA/concatallWCAfiltered2.vcf.gz $res/locator/concatWCAUlowread.vcf.gz  > $res/locator/concatallUhighreadR.vcf.gz

vcf-merge $res/locator/concat.allWCA/concatallWCAfiltered2.vcf.gz $res/locator/concatWCAUlowread.vcf.gz  > $res/locator/concatallUhighreadR.vcf.gz
exit

res=/work/idumville/pangolins/RAD/results
gzip $res/stacks*/populations*/*.vcf

exit