#!/bin/bash
#SBATCH --array=144-302%75
#SBATCH --error=/work/idumville/pangolins/RAD/bin/batch_outputs/%x%J%A_%a.error
#SBATCH --out=/work/idumville/pangolins/RAD/bin/batch_outputs/%x%J%A_%a.out
#SBATCH --job-name minmac2loc

#nofiltlocandPAs $nofiltlocator minmac3PAs  minmac2loc ALlocPAsminmac3us2
#1-365%50 0-5%5 6-364%75 2-364%362 2-364%25 1-40%40  1-365%50 0-32%32 0-365%5
data=$1
bin=$2
base_loc_dir=$3
res=/work/idumville/pangolins/RAD/results

echo module loading
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/bcftools-1.3.1
module load compiler/gcc-7.2.0
module load bioinfo/stacks-2.5

loc=concat
prefix=minmac2 #ALminmac3 #2 3

readarray -t VARS < $data/WCAurban.txt # $data/urbansamples2.txt #missingvcf.txt #urbansamples.txt #$data/WCAurban.txt #urbansamples2 is the ones that failed the first time due to oom
VAR=${VARS[$SLURM_ARRAY_TASK_ID]}
export VAR
echo $VAR

vcffile="/work/idumville/pangolins/RAD/results/locator/minmac2wRsalllineagesURpopmap_Doualacomb.vcf.gz"
#vcffile="/work/idumville/pangolins/RAD/results/locator/concatfull.vcf.gz"

oldmetadata="/work/idumville/pangolins/RAD/results/locator/metadata.txt"
rurlist="/work/idumville/pangolins/RAD/data/WCArurallist.txt" # "/work/idumville/pangolins/RAD/results/locatorAL/ruralsamples.txt" # 
urblist="/work/idumville/pangolins/RAD/data/WCAurban.txt" #"/work/idumville/pangolins/RAD/results/locatorAL/urbansamples.txt" # 

indv_loc_dir=$base_loc_dir/${VAR}_${prefix}_lociindv ; mkdir -p $indv_loc_dir

y=$(ls $indv_loc_dir/*vcf*)
p=$(ls $indv_loc_dir/*_predlocs.txt | wc -l)

if [ $p != 200 ]; then
  if [ -z "$y" ]
    then
    echo vcf production
    $bin/indvlocvcfs.sh $vcffile $indv_loc_dir $VAR $oldmetadata $rurlist
    rm $indv_loc_dir/locilist.txt #big file removal
    else
    echo no vcf production required
   fi 
  echo vcf is made
  echo locator
  $bin/iterateseedlocator.sh $VAR $prefix /work/idumville/pangolins/RAD/results $loc $indv_loc_dir/${VAR}lociindv.vcf.gz $indv_loc_dir/metadata.txt $indv_loc_dir
  echo removing the vcf
  rm $indv_loc_dir/${VAR}lociindv.vcf.gz
else
  echo this ${VAR} ${prefix} is done
fi

exit

#now we run populations
nt=4

loc=concat
stdir=$res/stacksautosomes
popdir=$stdir/populations${loc}.${VAR}_${prefix}_2 ; mkdir -p $popdir
popmap=$popdir/${VAR}stackspopmap.txt
awk '{print $1, "Train"}' $data/ruralsamples.txt | sed 's/ /\t/g' > $popmap
echo $VAR Test | sed 's/ /\t/g' >> $popmap
#popmap=$res/locator/${VAR}.allWCA00_lociindv/${VAR}stackspopmap.txt

echo running the populations
populations -V $indv_loc_dir/${VAR}lociindv.vcf.gz -M $popmap -O $popdir -t $nt --no-hap-exports

echo removing or zipping pop files
#rm ${popdir}/*log
rm ${popdir}/*distribs
gzip -f ${popdir}/*sumstats.tsv
rm ${popdir}/*hapstats.tsv
rm ${popdir}/*haplotypes.tsv
#rm $indv_loc_dir/${VAR}lociindv.vcf.gz



