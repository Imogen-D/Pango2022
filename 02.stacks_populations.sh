#!/bin/sh
nt=1 ; memo=1
nt=4 ; memo=10 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i


user=idumville
#user=jsalmona
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data

oldres=/work/jsalmona/pangolins/RAD/results
cov_dir=$res/coverage ; mkdir -p $cov_dir
stdir=$res/stacks ; mkdir -p $stdir
gdir=/work/jsalmona/pangolins/ref_genomes
genome=$gdir/Jaziri_pseudohap2_scaffolds_HiC


all_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-all
auto_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-autosomes
mt_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-mito
x_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-xchrom


mkdir -p $bin/batch_outputs

cd $bin
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

##make popmap file
#to be edited with population info later
#sed '1d' IF COPIED HEADERS
awk '{print $25, $26}' /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | sed 's/ /\t/' > $data/all.popmap

##################################
###### run gstacks ###############
##################################

nt=24 ; memo=98
loc=mito
stdir=$res/stacks${loc} ; mkdir -p $stdir
#popmap=$data/miseq_2021.popmap ; suffix='.rad.ptri.sorted_filtered.bam'

popmap=$data/all.popmap ; suffix=.rad.ptri.sorted_filtered.bam #sorted_filtered_${loc}.bam #'.rad.ptri.sorted_filtered.bam'

#sbatch --cpus-per-task=${nt} --mem=${memo}G  -J gs.miseq2021 -o $bin/batch_outputs/gs.miseq2021 $bin/run.gstacks.sh $nt $sb_dir $popmap $stdir $suffix
sbatch --cpus-per-task=${nt} --mem=${memo}G  -J gstacks.${loc}2021 -o $bin/batch_outputs/gstacks.${loc} $bin/run.gstacks.sh $nt $mt_sb_dir $popmap $stdir $suffix
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

##################################
###### run populations ###########
##################################
nt=1 ; memo=200
#popmap=$data/miseq_2021.popmap
loc=autosomes
prefix=alllineagesURpopmap_Doualacomb #cameroonURbiokoallreads #alllineagesUR #cameroonR #
popmap=$res/stacksautosomes/${prefix}popmap.txt #$data/all.popmap ; 
suffix=.sorted_filtered_${loc}.bam
stdir=$res/stacks${loc}
var=minmac2wRs #wRs04 #no_param #no_param  ##minmac2 #no_param #remember to change run.population.sh script to suit
popdir=$stdir/populations${loc}.${var}${prefix}vcf ; mkdir -p $popdir

r=0.4 ; R=0.4
sbatch --cpus-per-task=${nt} --mem=${memo}G  -J vcfmake.pop.${var}.${prefix} -o $bin/batch_outputs/gpop.${loc}.${var}${prefix}vcf $bin/run.population.sh $nt $popmap $stdir $popdir $r $R
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

#########rerunnning for lowread samples for Locator#########
prefix=WCAall #WCAallUhighR
loc=xchrom
suffix=.sorted_filtered_${loc}.bam
stdir=$res/stacks${loc}
var=mac2
popdir=$stdir/populations${loc}.${var}${prefix} ; mkdir -p $popdir
touch $popdir/${prefix}popmap.txt
popmap=$popdir/${prefix}popmap.txt ; rm $popmap

while read indv ; do
grep $indv /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | awk '{print $25, $39}'  | sed 's/ /\t/' >> $popmap
done < $data/WCAsamples.txt

nt=1 ; memo=200 #nt=30 ; memo=50


#r=0.4 ; R=0.4

sbatch --cpus-per-task=${nt} --mem=${memo}G  -J pop.${var}.${prefix} -o $bin/batch_outputs/gpop.${loc}.${var}${prefix} $bin/run.population.sh $nt $popmap $stdir $popdir # $r $R
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


######### rerunnning for indivudal samples #########
loc=autosomes #xchrom
suffix=.sorted_filtered_${loc}.bam
stdir=$res/stacks${loc}
prefix=allWCA

nt=1 ; memo=15     # X memo=10

r=0.4 ; R=0.4

rurlist=$data/rurallist.txt


while read indv ; do
popdir=$stdir/populations${loc}.${indv}wRs04 ; mkdir -p $popdir
indv_loc_dir=$res/locator/${indv}.${prefix}00_lociindv
rm $indv_loc_dir/${indv}stackspopmap.txt ; touch $indv_loc_dir/${indv}stackspopmap.txt
popmap=$indv_loc_dir/${indv}stackspopmap.txt
awk '{print $1}' $rurlist | while read rurindv; do 
printf "$rurindv Train\n" | sed 's/ /\t/g'  >> $indv_loc_dir/${indv}stackspopmap.txt
done
printf "$indv Test" | sed 's/ /\t/g'  >> $indv_loc_dir/${indv}stackspopmap.txt
sbatch --cpus-per-task=${nt} --mem=${memo}G  -J pop.${indv}.wRs04${loc} -o $bin/batch_outputs/gpop.${loc}.${indv}.wRs04 $bin/run.population.sh $nt $popmap $stdir $popdir $r $R
done < $data/WCAurban.txt



squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

loc=xchrom
stdir=$res/stacks${loc}
while read indv ; do
popdir=$stdir/populations${loc}.${indv}wRs04
gzip $popdir/*vcf
done < $data/WCAurban.txt

#non iterative but may work
nt=1 ; memo=25
sbatch --cpus-per-task=${nt} --mem=${memo}G -o $bin/batch_outputs/gpop.arrayminmac2 $bin/run.indv.popn.sh $res $data
## on Y200_1 + CAM004_1

nt=1 ; memo=25
sbatch --cpus-per-task=${nt} --mem=${memo}G -o $bin/batch_outputs/gpop.arrayseed $bin/run.indvseedpopn.sh $res $data

###removing when it failed
loc=autosomes
stdir=$res/stacks${loc}
while read indv ; do
popdir=$stdir/populations${loc}.${indv}minmac2wRs04 
#ls $popdir
tail $popdir/populations.log
rm -r $popdir
done < $data/WCAurban.txt


########## Extracting PA's #######
while read indv ; do
 echo  $indv all "$(grep -m 3 Test  /save/idumville/pangolins/RAD/results/stacksautosomes/populationsconcat.${indv}concat00seed/concatallWCAfiltered.p.log | tail -n 1 | awk '{print $10, $13}' | sed 's|/|\t|g' | sed 's/ /\t/g' | sed 's/;//')"  >> $res/privateallele.txt
  echo  $indv 90 "$(grep -m 3 Test /save/idumville/pangolins/RAD/results/stacksautosomes/populationsconcat.${indv}concat90seed/concatallWCA90filtered.p.log | tail -n 1 | awk '{print $10, $13}' | sed 's|/|\t|g' | sed 's/ /\t/g' | sed 's/;//')"  >> $res/privateallele.txt
  echo  $indv minmac3 "$(grep -m 3 Test /save/idumville/pangolins/RAD/results/stacksautosomes/populationsconcat.${indv}_minmac3*/*log | tail -n 1 | awk '{print $10, $13}' | sed 's|/|\t|g' | sed 's/ /\t/g' | sed 's/;//')"  >> $res/privateallele.txt
done < $data/WCAurban.txt
#number sites/#variant sites #polysites #PAs


sed -i '1 i\sample_type\tnew_name\tsites\tvariant\tpolymorphic\tPAs' $res/privateallele.txt


##################################
#### estimate error rate #########
##################################
#######haven't ran this yet

# recode plink
#module load bioinfo/plink-v1.90b5.3
#plink --allow-extra-chr --make-bed --out $stdir/populations/populations.plink.rcA --recodeA --file $stdir/populations/populations.plink
# filter vcf

loc=autosome #xchrom #autosome
var=minmac2 #no_param #minmac2
popdir=$stdir/populations${loc}.${var}

module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
# based on minimum read depth per genotype
for DP in $(seq 1 15) ; do
  vcftools --vcf $stdir/populations/populations.snps.vcf  --minDP $DP --recode --out $stdir/populations/populations.snps.dp.$DP
done
# based on minimum genotype quality score per genotype
for GQ in $(seq 5 5 40) ; do
  vcftools --vcf $stdir/populations/populations.snps.vcf  --minGQ $GQ --recode --out $stdir/populations/populations.snps.gq.$GQ
done
# run error rate estimation
module load system/R-3.5.2 libraries/gdal-2.3.0 libraries/proj-4.9.3 system/R-3.5.2 ; 
Rscript $bin/error_rates.r
