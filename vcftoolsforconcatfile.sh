#!/bin/sh

res=$1 ; echo $1
bin=$2 ; echo $2

echo loading modules
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/bcftools-1.3.1

prefix=minmac2wRsalllineagesURpopmap_Doualacomb
echo bgzipping
#bgzip -f -c $res/stacksautosomes/populationsautosomes.wRs04WCAallUhighR/populations.snps.vcf > $res/locator/autopopulations.snps.vcf.gz
#bgzip -f -c $res/stacksxchrom/populationsxchrom.wRs04WCAallUhighR/populations.snps.vcf > $res/locator/xchrompopulations.snps.vcf.gz

bgzip -f -c $res/stacksautosomes/populationsautosomes.${prefix}vcf/populations.snps.vcf > $res/locator/autosomes${prefix}populations.snps.vcf.gz
bgzip -f -c $res/stacksxchrom/populationsxchrom.${prefix}vcf/populations.snps.vcf > $res/locator/xchrom${prefix}populations.snps.vcf.gz


echo sorting
#vcf-sort $res/locator/autopopulations.snps.vcf.gz > $res/locator/autopopulationssort.snps.vcf
#vcf-sort $res/locator/xchrompopulations.snps.vcf.gz > $res/locator/xchrompopulationssort.snps.vcf
vcf-sort $res/locator/autosomes${prefix}populations.snps.vcf.gz | bgzip -f -c > $res/locator/auto${prefix}populationssort.snps.vcf.gz
vcf-sort $res/locator/xchrom${prefix}populations.snps.vcf.gz | bgzip -f -c > $res/locator/xchrom${prefix}populationssort.snps.vcf.gz

echo rezipping
#bgzip -f -c  $res/locator/autopopulationssort.snps.vcf > $res/locator/autopopulationssort.snps.vcf.gz
#bgzip -f -c $res/locator/xchrompopulationssort.snps.vcf > $res/locator/xchrompopulationssort.snps.vcf.gz


echo bcf indexing 
#bcftools index -t -f $res/locator/autopopulationssort.snps.vcf.gz
#bcftools index -t -f  $res/locator/xchrompopulationssort.snps.vcf.gz

bcftools index -t -f $res/locator/auto${prefix}populationssort.snps.vcf.gz
bcftools index -t -f $res/locator/xchrom${prefix}populationssort.snps.vcf.gz

echo tabix indexing
#tabix -p vcf $res/locator/autopopulationssort.snps.vcf.gz
#tabix -p vcf  $res/locator/xchrompopulationssort.snps.vcf.gz

tabix -p vcf $res/locator/auto${prefix}populationssort.snps.vcf.gz
tabix -p vcf $res/locator/xchrom${prefix}populationssort.snps.vcf.gz

echo concatenating
#bcftools concat $res/locator/xchrompopulationssort.snps.vcf.gz $res/locator/autopopulationssort.snps.vcf.gz > $res/locator/concatallUhighR.vcf.gz
bcftools concat $res/locator/xchrom${prefix}populationssort.snps.vcf.gz $res/locator/auto${prefix}populationssort.snps.vcf.gz > $res/locator/${prefix}.vcf

bgzip -f -c $res/locator/${prefix}.vcf  > $res/locator/${prefix}.vcf.gz 

vcftools --vcf $res/locator/${prefix}.vcf --plink --out $res/locator/${prefix}plink.plink

exit

echo unzipping and bgzipping
#gunzip -c $res/stacksautosomes/populationsautosomes.wRs04WCAallUhighR/populations.snps.vcf.gz | bgzip -f -c > $res/locator/autopopulations.snps.vcf.gz
#gunzip -c $res/stacksxchrom/populationsxchrom.wRs04WCAallUhighR/populations.snps.vcf.gz | bgzip -f -c > $res/locator/xchrompopulations.snps.vcf.gz

bgzip -f -c $res/stacksautosomes/populationsautosomes.wRs04WCAallUhighR/populations.snps.vcf > $res/locator/autopopulations.snps.vcf.gz
bgzip -f -c $res/stacksxchrom/populationsxchrom.wRs04WCAallUhighR/populations.snps.vcf > $res/locator/xchrompopulations.snps.vcf.gz

echo indexing
bcftools index -t -f $res/locator/autopopulations.snps.vcf.gz
bcftools index -t -f $res/locator/xchrompopulations.snps.vcf.gz

echo concatenating
bcftools concat $res/locator/xchrompopulations.snps.vcf.gz $res/locator/autopopulations.snps.vcf.gz > $res/locator/concatallUhighR.vcf.gz

exit

echo removing intemediates
rm $res/locator/autopopulations.snps.vcf.gz
rm $res/locator/xchrompopulations.snps.vcf.gz


exit 

#for 0% 95% filtering
loc=concat
prefix=allWCAlowreadU

base_loc_dir=$res/locator
vcffile=$res/locator/concatfull.vcf.gz #$res/stacks${loc}/populations${loc}.wRs04/populations.snps.vcf.gz 
nt=1 ; memo=6

VARS=("" 95) #"" 95
for mnum in "${VARS[@]}"; do
if [[ ! -z $mnum ]] ; then
mm=mm${mnum}
#if mm is a number, this filters by parameters in file (currently --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.75 --mac 3) then uses {prefix}.txt to extract samples by name 
sbatch --cpus-per-task=${nt} --mem=${memo}G  -J vcffile${loc}${prefix}${mm} -o $bin/batch_outputs/vcffile${loc}${prefix}${mm} $bin/vcftools.sh $vcffile $base_loc_dir $loc ${prefix} $base_loc_dir $mnum
else
mm=
#if mm is a number, this filters by parameters in file (currently --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.75 --mac 3) then uses {prefix}.txt to extract samples by name 
sbatch --cpus-per-task=${nt} --mem=${memo}G  -J vcffile${loc}${prefix}${mm} -o $bin/batch_outputs/vcffile${loc}${prefix}${mm} $bin/vcftools.sh $vcffile $base_loc_dir $loc ${prefix} $base_loc_dir $mnum
fi
done

exit


