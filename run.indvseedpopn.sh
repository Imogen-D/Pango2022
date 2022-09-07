#!/bin/bash
#SBATCH --job-name indvpopnsseed
#SBATCH --array=1-365%20

##give sbatch these
res=$1
nt=1
data=$2
r=0.4

module load compiler/gcc-7.2.0
module load bioinfo/stacks-2.5

readarray -t VARS < $data/WCAurban.txt
VAR=${VARS[$SLURM_ARRAY_TASK_ID]}
export VAR

prefix=allWCA
loc=concat
stdir=$res/stacksautosomes

indv_loc_dir=$res/locator/${VAR}.${prefix}00_lociindv
base_loc_dir=$res/locator
vcffile=$base_loc_dir/concat.allWCAmm90/concatallWCA90filtered.vcf.gz

popdir=$stdir/populations${loc}.${VAR}concat90seed ; mkdir -p $popdir
popmap=$indv_loc_dir/${VAR}stackspopmap.txt
populations -V $vcffile -M $popmap -O $popdir -t $nt
#rm ${popdir}/*log
rm ${popdir}/*distribs
#rm ${popdir}/*sumstats.tsv
rm ${popdir}/*hapstats.tsv
rm ${popdir}/*haplotypes.tsv


vcffile=$base_loc_dir/concat.allWCA/concatallWCAfiltered.vcf.gz

popdir=$stdir/populations${loc}.${VAR}concat00seed ; mkdir -p $popdir
popmap=$indv_loc_dir/${VAR}stackspopmap.txt
populations -V $vcffile -M $popmap -O $popdir -t $nt
#rm ${popdir}/*log
rm ${popdir}/*distribs
#rm ${popdir}/*sumstats.tsv
rm ${popdir}/*hapstats.tsv
rm ${popdir}/*haplotypes.tsv


