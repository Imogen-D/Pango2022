#!/bin/bash
#SBATCH --job-name indvpopnsMINMAC2
#SBATCH --array=1-20%20
#SBATCH --error=/work/idumville/pangolins/RAD/bin/batch_outputs/%x%J%A_%a.error
#SBATCH --out=/work/idumville/pangolins/RAD/bin/batch_outputs/%x%J%A_%a.out

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

loc=concat
stdir=$res/stacksautosomes
base_loc_dir=$res/locator
vcffile=$base_loc_dir/${VAR}_minmac2_lociindv/${VAR}lociindv.vcf.gz
popdir=$stdir/populations${loc}.${VAR}_minmac2 ; mkdir -p $popdir
prefix=allWCA; indv_loc_dir=$res/locator/${VAR}.${prefix}00_lociindv ; popmap=$indv_loc_dir/${VAR}stackspopmap.txt

populations -V $vcffile -M $popmap -O $popdir -t $nt --no-hap-exports
#rm ${popdir}/*log
rm ${popdir}/*distribs
gzip ${popdir}/*sumstats.tsv
rm ${popdir}/*hapstats.tsv
rm ${popdir}/*haplotypes.tsv


exit

##OLD WAY WITH STACKS DIRECTORY
#!/bin/bash
#SBATCH --job-name indvpopnsminmac2
#SBATCH --array=1-2%2

#1-365%75

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


loc=autosomes #xchrom
suffix=.sorted_filtered_${loc}.bam
stdir=$res/stacks${loc}
prefix=allWCA
rurlist=$data/WCArurallist.txt

popdir=$stdir/populations${loc}.${VAR}minmac2wRs04 ; mkdir -p $popdir
indv_loc_dir=$res/locator/${VAR}.${prefix}00_lociindv

rm $indv_loc_dir/${VAR}stackspopmap.txt ; touch $indv_loc_dir/${VAR}stackspopmap.txt
popmap=$indv_loc_dir/${VAR}stackspopmap.txt
awk '{print $1}' $rurlist | while read rurindv; do 
printf "$rurindv Train\n" | sed 's/ /\t/g'  >> $indv_loc_dir/${VAR}stackspopmap.txt
done
printf "$VAR Test" | sed 's/ /\t/g'  >> $indv_loc_dir/${VAR}stackspopmap.txt

populations -P $stdir/ -M $popmap -O $popdir -t $nt --min-mac 2 -r $r #-p 1 
#rm ${popdir}/*log
rm ${popdir}/*distribs
gzip ${popdir}/*sumstats.tsv
rm ${popdir}/*hapstats.tsv
rm ${popdir}/*haplotypes.tsv

loc=xchrom #xchrom
suffix=.sorted_filtered_${loc}.bam
stdir=$res/stacks${loc}
prefix=allWCA
rurlist=$data/rurallist.txt

popdir=$stdir/populations${loc}.${VAR}minmac2wRs04 ; mkdir -p $popdir
indv_loc_dir=$res/locator/${VAR}.${prefix}00_lociindv

#rm $indv_loc_dir/${VAR}stackspopmap.txt ; touch $indv_loc_dir/${VAR}stackspopmap.txt
popmap=$indv_loc_dir/${VAR}stackspopmap.txt
#awk '{print $1}' $rurlist | while read rurindv; do 
#printf "$rurindv Train\n" | sed 's/ /\t/g'  >> $indv_loc_dir/${VAR}stackspopmap.txt
#done
#printf "$VAR Test" | sed 's/ /\t/g'  >> $indv_loc_dir/${VAR}stackspopmap.txt

populations -P $stdir/ -M $popmap -O $popdir -t $nt --min-mac 2 -r $r #-p 1 
#rm ${popdir}/*log
rm ${popdir}/*distribs
gzip ${popdir}/*sumstats.tsv
rm ${popdir}/*hapstats.tsv
rm ${popdir}/*haplotypes.tsv