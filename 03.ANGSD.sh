#!/bin/sh
#SBATCH -p workq
nt=4 ; memo=10 
nt=1 ; memo=4 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i #--x11 

# define paths

user=idumville
#user=jsalmona
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data

oldres=/work/jsalmona/pangolins/RAD/results
olddata=/work/jsalmona/pangolins/RAD/data

all_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-all
auto_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-autosomes
mt_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-mito
x_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-xchrom

cov_dir=$res/coverage ; mkdir -p $cov_dir
stdir=$res/stacks ; mkdir -p $stdir

gdir=/work/jsalmona/pangolins/ref_genomes
genome=$gdir/Jaziri_pseudohap2_scaffolds_HiC

angsd_dir=$res/angsd ; mkdir -p $angsd_dir
ngsadmix_dir=$res/ngsadmix ; mkdir -p $ngsadmix_dir
cd $bin

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

# grab sample list and number of samples
awk '{print $25}' /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt > $data/pango_sample_list.txt
N=$(cat $data/pango_sample_list.txt | wc -l | cut -d " " -f1)
echo number of individuals = $N

############STILL TO DO? #############
###### estimate and plot min and max coverage per individual and globally ######
  loc=autosomes
  sample_file=pango.rad.all_$loc
  echo plot loci coverage distribution and estimate min and max thresholds: min=N*4 max=quantile0.9975
  rm $cov_dir/minmax_${prefix}.txt ; touch $cov_dir/minmax_${prefix}.txt
  # module load system/R-3.4.3 ; Rscript --vanilla $BIN/indv_cov_plot.r $nt $sample_file $N $prefix > cov.$prefix
  nt=2
  memo=4
  sbatch --cpus-per-task=${nt} --mem=${memo}G  -J cov.$prefix. -o $BIN/cov.$prefix $BIN/cov_R.sh $nt $sample_file $N $prefix $BIN $memo
squeue -u $user 

#######################################################################################################
####### estimate genotype likelihood with angsd and store them in a beagle file for ngsadmix ########## 
#######################################################################################################
##rerunning with 
#manually removed CA, Gha, WAf, L27, Y360, Y199, now rmeoving cov above 3
awk '$5>=3' cleanedsamples1.txt | awk '{print $1}' > $res/cleanedsamples.txt ; rm $res/cleanedsamples1.txt
#from now on put in SEPERATE angsd folders? Overwritten old angsd so what currently in /angsd/ is just autosomes
prefix=rmCADGWaf #parameters used to be prefix=rad.ptri
loc=xchrom
angsd_dir=$res/angsd.${prefix}.${loc} ; mkdir -p $angsd_dir



nt=30 ; memo=100

sample_file=pango.rad.${loc}_$prefix


genome=$gdir/Jaziri_pseudohap2_scaffolds_HiC.fasta #KP306514.1.fasta #Jaziri_pseudohap2_scaffolds_HiC.fasta
  rm $res/${sample_file}.bamlist # create bam list #$prefix.
  rm $cov_dir/ind_cov.${sample_file}.txt #.${prefix}
  cat $res/cleanedsamples.txt | while read indv ; do #$data/pango_sample_list.txt
   #only for xchrom remove 199_1
   # if [[ "$indv" == 'Y199_1' ]]; then
    #  continue
    #fi
    #CHANGE SB DIR
    echo $x_sb_dir/${indv}.sorted_filtered_${loc}.bam >> $res/${sample_file}.bamlist #.rad.ptri.sorted_filtered.bam >> $res/${sample_file}.bamlist
    awk -v var="$indv" '$1 == var' $cov_dir/RAD_${loc}_cov_stats.txt >> $cov_dir/ind_cov.${sample_file}.txt ; done 
  head $res/${sample_file}.bamlist
  head $cov_dir/ind_cov.${sample_file}.txt
  
  #for XCHROM Y199_1 MANUALLY REMOVED to check later

#this stat doesn't seem right so not using this parameter for now
#  maxdepthind=$(cat $cov_dir/ind_cov.${sample_file}_${loc}.txt | cut -f5 | sort -n | tail -1) ; echo "maxdepthind = " $maxdepthind
#maxdepthin=500
  #this used to divide by 4 and had opp gmin gmax but unsure why? Now just as is 
  gmax=$(cat $cov_dir/ind_cov.${sample_file}.txt | cut -f5 | paste -sd+ | bc) ; echo "gmax = " $gmax
  
  #gmin is essentially just two times no of samples
  gmin=$(cat $cov_dir/ind_cov.${sample_file}.txt | cut -f4 | paste -sd+ | bc) ; echo "gmin = " $gmin
  minind=1 ; echo "minind = " $minind
  mindepthind=1 #for initial run
  cat $res/cleanedsamples.txt | wc -l #cat $data/pango_sample_list.txt | wc -l
  cat $res/${sample_file}.bamlist | wc -l 
  cat $cov_dir/ind_cov.${sample_file}.txt | wc -l  #${prefix}.
 
#for mito also changed angsd.sh to not use gmax TO CHANGE BACK


#for mit
#gmax=20000
gmin=1
minind=1
mindepthind=1
maxdepthind=1 #not used
sample_file=pango.rad.${loc}_$prefix
bamlist=$res/${sample_file}.bamlist
#BE CAREFUL ABOUT SB_DIR BECAUSE IF TO FIND BAMLIST IS NOW IN IDUMVILLE RES ; IF TO FIND BAMS IS IN AUTO_SB_DIR
  sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ang.$prefix.$loc -o $bin/batch_outputs/ang-$prefix-$loc $bin/ANGSD.sh $nt $bamlist $genome $angsd_dir $prefix $sample_file $gmin $gmax $maxdepthind $mindepthind $minind
  #$bin/angsd.mms.sh $nt $sb_dir $genome $angsd_dir $prefix $sample_file $gmin $gmax $maxdepthind $mindepthind $minind
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"
#######################################################################################################
  
########
prefix=rmCADGWaf #parameters used to be prefix=rad.ptri
loc=autosomes
sample_file=pango.rad.${loc}_$prefix
angsd_dir=$res/angsd.${prefix}.${loc} ; mkdir -p $angsd_dir
ls -lh $angsd_dir 
# how many markers are included 
zcat $angsd_dir/$prefix.${sample_file}.beagle.gz | wc -l &
zcat $angsd_dir/$prefix.${sample_file}.beagle.gz | wc -l  > $angsd_dir/$prefix.${sample_file}.angsd.nb.of.loci #store this number &

#######################################################################################################
####### run ngsadmix                                                                         ########## 
#######################################################################################################
prefix=rmCADGWaf #parameters used to be prefix=rad.ptri
loc=xchrom
sample_file=pango.rad.${loc}_$prefix
angsd_dir=$res/angsd.${prefix}.${loc} ; mkdir -p $angsd_dir

ngsadmix_dir=$res/ngsadmix.$loc.${prefix} ; mkdir -p $ngsadmix_dir

#  maxdepthind=$(cat $cov_dir/ind_cov.${sample_file}_${loc}.txt | cut -f5 | sort -n | tail -1) ; echo "maxdepthind = " $maxdepthind
maxdepthind=1000 #placeholder value
  gmax=$(cat $cov_dir/ind_cov.${sample_file}.txt | cut -f5 | paste -sd+ | bc) ; echo "gmax = " $gmax
  
  #gmin is essentially just two times no of samples
  gmin=$(cat $cov_dir/ind_cov.${sample_file}.txt | cut -f4 | paste -sd+ | bc) ; echo "gmin = " $gmin
  #gmax=1000 #not used; no indcov file for mito
  #gmin=1000 #not used; no indcov file for mito
  minind=1 ; echo "minind = " $minind
  mindepthind=1 #for initial run
 
  for K in $(seq 1 10) ; do for seed in $(seq 1 20) ; do 
    nt=16 ; memo=32 ; if [[ "$K" == "1" ]] ; then nt=8 ; memo=20 ; fi #doubled these for autosomes was 4;8 and 1;2 ish     
        echo $nt $prefix $K $seed $sample_file $gmin $gmax $minind $maxdepthind $bin $ngsadmix_dir $angsd_dir ;  
        sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ngad.$sample_file.${K}.$seed -o $bin/batch_outputs/ngad.$sample_file.${K}.$seed.2 $bin/ngsadmix.sh $nt $prefix $K $seed $sample_file $gmin $gmax $minind $maxdepthind $bin $ngsadmix_dir $angsd_dir;
done ; done 
# monitor your jobs :
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"
# monitor your ngsadmix output files
ls -lh $ngsadmix_dir/*/*.qopt | wc -l

#######################################################################################################
####### cleaning beagle file to rerun ngsadmix                                          ########## 
#######################################################################################################
prefix=rad.ptri.pango.rad.all_
loc=xchrom

ngsadmix_dir=/work/idumville/pangolins/RAD/results/ngsadmix.${loc}.lowreads/ ; mkdir -p $ngsadmix_dir

nt=5 ; memo=10

output=allhighread_wobio #ruralhighread_wobio #list to keep #ruralhighread_wobiogab #ruralhighread
input=$res/ngsadmix.autosomes/lowurbanreads_wobio.txt #list to remove
sbatch --cpus-per-task=${nt} --mem=${memo}G  -J cutting_beagle -o $bin/batch_outputs/beaglefilt${loc}${output} $bin/cuttingbeagle.sh ${loc} $res $ngsadmix_dir $output $input

#for IBD files
rm $res/ngsadmix.autosomes.lowreads/concatruralhighread_wobiobeagle.gz
cp $res/ngsadmix.autosomes.lowreads/autosomesruralhighread_wobiobeagle.gz $res/ngsadmix.autosomes.lowreads/concatruralhighread_wobiobeagle.gz

gunzip -c $res/ngsadmix.xchrom.lowreads/xchromruralhighread_wobiobeagle.gz | sed  '1d' | gzip -c >> $res/ngsadmix.autosomes.lowreads/concatruralhighread_wobiobeagle.gz

# monitor your jobs :
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


#######################################################################################################
####### new ngsadmix with high read samples                                         ########## 
#######################################################################################################
loc=xchrom

prefix=lowreads
ngsadmix_dir=$res/ngsadmix.$loc.$prefix ; mkdir -p $ngsadmix_dir
listsamples=$res/ngsadmix.${loc}.lowreads/${loc}ruralhighreadsamples.txt #$ngsadmix_dir
beagle=$ngsadmix_dir/${loc}ruralhighreadbeagle.gz

fileprefix=ruralhighread #parameters used to be prefix=rad.ptri
sample_file=pango.rad.${loc}
angsd_dir=$res/angsd.${prefix}.${loc} ; mkdir -p $angsd_dir

#new ngsadmix folder for these rural
ngsadmix_dir=$res/ngsadmix.$loc.$fileprefix ; mkdir -p $ngsadmix_dir

maxdepthind=1000 #placeholder value, unused
  gmax=$(cat $cov_dir/ind_cov.pango.rad.all_autosomes.txt | cut -f5 | paste -sd+ | bc) ; echo "gmax = " $gmax
   #gmin is essentially just two times no of samples
  gmin=$(cat ${listsamples} | awk 'END {print NR*2}') ; echo "gmin = " $gmin
  minind=1 ; echo "minind = " $minind
  mindepthind=1 #for initial run
  
  for K in $(seq 1 10) ; do for seed in $(seq 1 20) ; do 
    nt=4 ; memo=10 ; if [[ "$K" == "1" ]] ; then nt=1 ; memo=2 ; fi #doubled these for autosomes to 16, 32 was 4;8 and 1;2 ish     
        echo $nt $fileprefix $K $seed $sample_file $gmin $gmax $minind $maxdepthind $bin $ngsadmix_dir $angsd_dir ;  
        sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ngad.$sample_file.${K}.$seed -o $bin/batch_outputs/ngad.${fileprefix}${sample_file}.${K}.$seed $bin/ngsadmix_newbeagle.sh $nt $fileprefix $K $seed $sample_file $gmin $gmax $minind $maxdepthind $bin $ngsadmix_dir $angsd_dir $beagle;
done ; done 

# monitor your jobs :
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"
# monitor your ngsadmix output files
ls -lh $ngsadmix_dir/*/*.qopt | wc -l


#######################################################################################################
###### produce likelihood value file ######

loc=autosomes
prefix=rad.ptri
sample_file=pango.rad.all_${loc}

ngsadmix_dir=$res/ngsadmix.${loc}
rm $ngsadmix_dir/$prefix.${sample_file}_likevalues.txt ; touch $ngsadmix_dir/$prefix.${sample_file}_likevalues.txt
for K in $(seq 1 10) ; do for seed in $(seq 1 20) ; do 
  grep best $ngsadmix_dir/$prefix.${sample_file}_k${K}/$prefix.${sample_file}_k${K}_s${seed}.log | awk '{print $K}' | cut -d'=' -f2- | sort -g | sed "s/after/$K/g" | sed "s/iterations/$seed/g" >> $ngsadmix_dir/$prefix.${sample_file}_likevalues.txt
done ; done ; 
cat $ngsadmix_dir/$prefix.${sample_file}_likevalues.txt

loc=xchrom
prefix=lowreads #prefix=rad.ptri
sample_file=pango.rad.${loc}
ngsadmix_dir=$res/ngsadmix.${loc}.$prefix
rm $ngsadmix_dir/$prefix.${sample_file}_likevalues.txt ; touch $ngsadmix_dir/$prefix.${sample_file}_likevalues.txt
for K in $(seq 1 20) ; do for seed in $(seq 1 20) ; do 
  grep best $ngsadmix_dir/$prefix.${sample_file}_k${K}/$prefix.${sample_file}_k${K}_s${seed}.log | awk '{print $K}' | cut -d'=' -f2- | sort -g | sed "s/after/$K/g" | sed "s/iterations/$seed/g" >> $ngsadmix_dir/$prefix.${sample_file}_likevalues.txt
done ; done ; 
cat $ngsadmix_dir/$prefix.${sample_file}_likevalues.txt

ngsadmix_dir=$res/ngsadmix.autosomes.lowreads
rm $ngsadmix_dir/$prefix.pango.rad.autosomes_likevalues.txt ; touch $ngsadmix_dir/$prefix.pango.rad.autosomes_likevalues.txt
for K in $(seq 1 10) ; do for seed in $(seq 1 20) ; do 
  grep best $ngsadmix_dir/$prefix.pango.rad.autosomes_k${K}/$prefix.pango.rad.autosomes_k${K}_s${seed}.log | awk '{print $K}' | cut -d'=' -f2- | sort -g | sed "s/after/$K/g" | sed "s/iterations/$seed/g" >> $ngsadmix_dir/$prefix.pango.rad.autosomes_likevalues.txt
done ; done ; 
cat $ngsadmix_dir/$prefix.${sample_file}_likevalues.txt

#######################################################################################################
########## plot and geoplot ngsadmix results ##############
covtable=$olddata/RAD_pango_recap.txt
popmap=$data/all.popmap
popmap2=$res/ngsadmix.${loc}/sampledata${loc}.txt ; rm $popmap2; touch $popmap2

#make new file for coverage with nloci/n_used_fw/mean_cov_ns/qmap
##for mito
loc=mito
covtable=$olddata/RAD_pango_recap.txt
popmap=$data/all.popmap
popmap2=$res/ngsadmix.${loc}/sampledata${loc}.txt ; rm $popmap2; touch $popmap2
awk '{print $1}' $popmap | while read indv; do
x=$(grep -w $indv $popmap)
y=$(grep -w $indv $covtable)
printf "$x $y\n" | awk '{print $1, $41, $39, $40, $38, $33, $8, $9, $30, $34}' | sed 's/ /\t/g' >> $popmap2 #$res/ngsadmix.${loc}/sampledata${loc}.txt
done
##for xchrom
loc=xchrom
covtable=$olddata/RAD_pango_recap.txt
popmap=$data/all.popmap
popmap2=$res/ngsadmix.${loc}/sampledata${loc}.txt ; rm $popmap2; touch $popmap2
#samples=$res/cleanedsamples.txt
awk '{print $1}' $popmap | while read indv; do #samples
x=$(grep -w $indv $popmap)
y=$(grep -w $indv $covtable)
printf "$x $y\n" | awk '{print $1, $41, $36, $37, $35, $31, $8, $9, $30, $34}' | sed 's/ /\t/g' | sed 's/,//g' >> $popmap2 #$res/ngsadmix.${loc}/sampledata${loc}.txt
done
##for autosomes
loc=autosomes
covtable=$olddata/RAD_pango_recap.txt
popmap=$data/all.popmap
popmap2=$res/ngsadmix.${loc}.$prefix/sampledata${loc}.txt ; rm $popmap2; touch $popmap2
samples=$res/ngsadmix.${loc}.$prefix/highreadsamples.txt
awk '{print $1}' $samples | while read indv; do #$popmap samples
x=$(grep -w $indv $popmap)
y=$(grep -w $indv $covtable)
printf "$x $y\n" | awk '{print $1, $41, $22, $23, $25, $32, $8, $9, $30, $34, $42}' | awk '{$12 = $11 / $5}1' | sed 's/ /\t/g' | sed 's/,//g' >> $popmap2 #$res/ngsadmix.${loc}/sampledata${loc}.txt 
done  



##for lowread xchrom
loc=xchrom
covtable=$olddata/RAD_pango_recap.txt
popmap=$data/all.popmap
popmap2=$res/ngsadmix.${loc}.$prefix/sampledata${loc}.txt ; rm $popmap2; touch $popmap2
samples=$res/ngsadmix.autosomes.$prefix/highreadsamples.txt
awk '{print $1}' $samples | while read indv; do #$popmap samples
x=$(grep -w $indv $popmap)
y=$(grep -w $indv $covtable)
printf "$x $y\n" | awk '{print $1, $41, $36, $37, $35, $31, $8, $9, $30, $34, $42}' | awk '{$12 = $11 / $5}1' | sed 's/ /\t/g' | sed 's/,//g' >> $popmap2 #$res/ngsadmix.${loc}/sampledata${loc}.txt  #40 for locality, 27 for popn
done  


##into geoplot script

#where is ID/locality/nloci/usedfwread/meancov/qmapreads/lat/long/lineage/market/He/Henormed
#could replace this script later with just using RAD_pango_recap not popmap but ensures consistency between pop of gstacks atm?

#covtable=$cov_dir/RAD_${loc}_cov_stats.txt


for loc in xchrom ; do #autosomes
  prefix=rad.ptri
  order_file=$data/pango_sample_list.txt
  sample_file=pango.rad.all_$loc
  ngsadmix_dir=$res/ngsadmix.$loc
  module load system/R-3.4.3
  Rscript --vanilla $bin/ngsadmix_plots.R $popmap2 $prefix.${sample_file} $ngsadmix_dir $order_file
  module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3
  Rscript --vanilla $bin/ngsadmix_geoplots.R $popmap2 $prefix.${sample_file} $ngsadmix_dir $order_file
done


#######################################################
###################### PCANGSD ########################
#######################################################

# for concatenated data, use loc = autosomes but change beagle and output name $res/ngsadmix.autosomes.lowreads/concatruralhighread_wobiobeagle.gz


### making new infofile
loc=autosomes #xchrom #autosomes
prefix=lowreads
fileprefix=ruralhighread_wobio #ruralhighread_wobio #gab


#if not run ngsadmix, manually make $res/ngsadmix.${loc}.$prefix/${loc}${fileprefix}samples.txt ie for allhighread_wobio

covtable=$olddata/RAD_pango_recap.txt
popmap=$data/all.popmap
popmap2=$res/ngsadmix.${loc}.$prefix/sampledata${loc}${fileprefix}.txt ; rm $popmap2; touch $popmap2
samples=$res/ngsadmix.${loc}.$prefix/${loc}${fileprefix}samples.txt #highreadsamples.txt #${loc}#${fileprefix}samples.txt
awk '{print $1}' $samples | while read indv; do #$popmap samples
x=$(grep -w $indv $popmap)
y=$(grep -w $indv $covtable)
printf "$x $y\n" | awk '{print $1, $41, $28}' | sed 's/ /\t/g' | sed 's/,//g' >> $popmap2 #$res/ngsadmix.${loc}/sampledata${loc}.txt  #indv, locality, popn
done  
#| awk '{print $1, $41, $28}'  

###Manually edited sampledataautosomesruralhighread_wobio.txt to include number for NGSADMIX colouring

sample_file=pango.rad.${loc} #pango.rad.all_$loc # with merged individual replicates

ngsadmix_dir=$res/ngsadmix.$loc.lowreads ; mkdir -p $ngsadmix_dir
#beaglefile=$ngsadmix_dir/${loc}${fileprefix}beagle.gz 
beaglefile=$ngsadmix_dir/concat${fileprefix}beagle.gz

#$ngsadmix_dir/lowreads_${loc}_beagle.gz #unfilteredallautosomes.beagle.gz # lowreads_${loc}_beagle.gz #$angsd_dir/$prefix.$sample_file.beagle.gz

#pcangsdir=$res/pcangsd.${loc}.${fileprefix} ; mkdir -p $pcangsdir
pcangsdir=$res/pcangsd.concat.${fileprefix} ; mkdir -p $pcangsdir

outfile=$pcangsdir/$fileprefix.$sample_file.all
infofile=$ngsadmix_dir/sampledata${loc}${fileprefix}.txt #res/ngsadmix.${loc}/sampledata${loc}.txt

nt=1 ; memo=4

for npc in all; do # 2 4 6 8  for making distance tables for pcangsd
outfile=$pcangsdir/$fileprefix.$sample_file.$npc
  echo $bin/batch_outputs/pc.$npc.$fileprefix.$sample_file $bin/pcangsd_pango.sh $beaglefile $outfile $nt $npc $pcangsdir $bin $infofile $loc
  sbatch --cpus-per-task=${nt} --mem=${memo}G  -J pc.$npc.$fileprefix.$sample_file -o $bin/batch_outputs/pc.$npc.$fileprefix.$sample_file.concat $bin/pcangsd_pango.sh $beaglefile $outfile $nt $npc $pcangsdir $bin $infofile $loc #sbatch --cpus-per-task=${nt} --mem=${memo}G  -J pc.$npc.$fileprefix.$sample_file -o $bin/batch_outputs/pc.$npc.$fileprefix.$sample_file
  #for ib in 1 2 3; do
   #$bin/pcangsd_inbreed_pango.sh $beaglefile $outfile.ib${ib} $nt $npc $pcangsdir $ib $bin $infofile $loc #sbatch --cpus-per-task=${nt} --mem=${memo}G  -J pc.$fileprefix.$sample_file -o $bin/pc.$fileprefix.$sample_file
#done 
done
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

# pcangsd + imbreeding + kinship


#######################################################
###################### IBD ############################
#######################################################
#making sample data frame with latlong


loc=autosomes #xchrom #autosomes
prefix=lowreads
fileprefix=ruralhighread_wobio #gab
IBD_dir=$res/IBD.${fileprefix} ; mkdir -p $IBD_dir

covtable=$olddata/RAD_pango_recap.txt
popmap=$data/all.popmap
popmap2=$IBD_dir/geodata${fileprefix}.txt ; rm $popmap2; touch $popmap2
samples=$res/ngsadmix.${loc}.$prefix/${loc}${fileprefix}samples.txt
awk '{print $1}' $samples | while read indv; do #$popmap samples
x=$(grep -w $indv $popmap)
y=$(grep -w $indv $covtable)
printf "$x $y\n" | awk '{print $1, $8, $9}'  | sed 's/ /\t/g' | sed 's/,//g' >> $popmap2 #$res/ngsadmix.${loc}/sampledata${loc}.txt  #indv, locality, popn
done  

DGen_concat=$res/pcangsd.concat.ruralhighread_wobio/ruralhighread_wobio.pango.rad.autosomes.allcov_euc_dist.rda
DGen_x=$res/pcangsd.xchrom.ruralhighread_wobio/ruralhighread_wobio.pango.rad.xchrom.allcov_euc_dist.rda
DGen_auto=$res/pcangsd.autosomes.ruralhighread_wobio/ruralhighread_wobio.pango.rad.autosomes.allcov_euc_dist.rda
Cam_samples=$res/IBD.ruralhighread_wobio/cameroon.txt
samplelist=$res/ngsadmix.autosomes.lowreads/autosomesruralhighread_wobiosamples.txt

#running Rscript
module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; Rscript $bin/IBD_plotting.R $IBD_dir geodata${fileprefix}.txt $samplelist $DGen_concat $DGen_x $DGen_auto $Cam_samples


###########################################################################
########## generate unfolded SFS for pogen He measures ##########
########## at the INDIVIDUAL level                               ##########
###########################################################################
prefix="rad.ptri"
loc=autosomes
hedir=$res/${loc}_indv_he ; mkdir -p $hedir


angsd_dir=$res/angsd ; mkdir -p $angsd_dir

#cat $data/pango_sample_list.txt |  while read subdata ; do
#cat $res/missinghe.txt |  while read subdata ; do
#subdata=DlaA30_3
    sample_file=indv
    rm $hedir/$prefix.$subdata.${sample_file}.bamlist
    ls $auto_sb_dir/${subdata}.sorted_filtered_${loc}.bam >> $hedir/$prefix.$subdata.${sample_file}.bamlist
    #ls $sb_dir/${subdata}_rmdups_noMTnoCP_${prefix}_sorted.bam > $sb_dir/$prefix.$subdata.${sample_file}.bamlist  # create bam list
    #cat $hedir/$prefix.$subdata.${sample_file}.bamlist
    #nt=4 ; memo=10 ; N=1 ; gmin=na ; gmax=na ; mindepthind=1 ; minind=1 ; genome=$denovoref
    nt=6 ; memo=20 ; N=1 ; gmin=na ; gmax=na ; mindepthind=1 ; minind=1 ; genome=$gdir/Jaziri_pseudohap2_scaffolds_HiC.fasta #2 4
    maxdepthind=$(cat $cov_dir/ind_cov.pango.rad.all_${loc}.txt | cut -f5 | sort -n | tail -1) ; echo "maxdepthind = " $maxdepthind
    echo $nt $auto_sb_dir $genome $angsd_dir $prefix $sample_file $gmin $gmax $maxdepthind $mindepthind $minind $subdata
#estimate genotype likelihood for sfs
    sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ang.he.$prefix.$subdata -o $bin/batch_outputs/ang.he.ind.$prefix.$subdata.$loc $bin/angsd.he.indv.sh $nt $hedir $genome $angsd_dir $prefix $sample_file $gmin $gmax $maxdepthind $mindepthind $minind $subdata $hedir/$prefix.$subdata.${sample_file}.bamlist
done
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R" 

loc=xchrom

#rm $res/${loc}hes.txt #$res/${loc}_indv_he/all_indv_he_${loc}.txt
cat $data/pango_sample_list.txt |  while read indv ; do

hedir=$res/${loc}_indv_he

he=$(cat $hedir/$prefix.$indv.indv.un.he)
echo $he
ho=$(cat $hedir/$prefix.$indv.indv.un.ho | sed '3q;d') # $res/ngsadmix.${loc}/sampledata${loc}.txt
x=$(grep -w $indv $res/ngsadmix.${loc}/sampledata${loc}.txt) #popmap2
echo $x

nsite=$(sed -n 1p $hedir/*${indv}.indv.un.ho)
nsites=$(sed -n 2p $hedir/*${indv}.indv.un.ho | awk -F"e" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}')
echo $nsite $nsites
VAR=$(echo "scale=10; $nsite/$nsites" | bc)

#printf "${x}\t${he}\t${VAR}\n" | sed 's/ /\t/g' >> $res/${loc}hes.txt # $res/${loc}_indv_he/all_indv_he_${loc}.txt  #\t${ho}
done



#now into plot_he.R
#where is ID/popn/nloci/usedfwread/meancov/qmapreads/lat/long/he/ho

############################################################################
############################################################################


########### Pop / site diversity estimates, He theta ... ########
# for localities
for subdata in $(cat $DATA/data.info | tail -n +2 | cut -f 3 | sort | uniq )  ; do echo $subdata
  sample_file="spini_all.$subdata" ; prefix="spini" ; genome=$denovoref;     nt=4 ; memo=8
    sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ang.th.$prefix.$subdata -o $BIN/temp_out_files/ang.th.$subdata $BIN/angsd_theta.spini.sh $nt $sb_dir $genome $angsd_dir $prefix $sample_file $subdata
done 
squeue -u $user


########### FST #############
# estimate fst among sampling sites
fst_dir=$RES/fst ; mkdir $fst_dir
nt=4
memo=8
prefix="spini"
boot=100 # these are actually permutations not bootstrap
#cat $DATA/data.info.2 | tail -n +2 | cut -f 3 | sort | uniq > $DATA/sampling.sites
#for pop1 in Benanofy ; do #$(cat $DATA/data.info.2 | tail -n +2 | cut -f 3 | sort | uniq | head -n -1)  ; do 
#  for pop2 in Solaniampilana ; do #$(cat $DATA/data.info.2 | tail -n +2 | cut -f 3 | sort | uniq | tail -n +1)  ; do 
cat $DATA/forest_pair_list.2 | tail -n +11 | while read pop1 pop2 ; do
  echo $pop1 $pop2
    if [[ "$pop1" == "$pop2" ]] ; then echo "pop1=pop2 moving to next pair" ; else
    if [[ "$pop1" == "Ampondrabe" ]] ; then echo "pop1=Ampondrabe moving to next pair" ; else
    if [[ "$pop2" == "Ampondrabe" ]] ; then echo "pop2=Ampondrabe moving to next pair" ; else
    #sbatch --cpus-per-task=${nt} --mem=${memo}G -J ang.Fst.$pop1.$pop2 -o $BIN/temp_out_files/ang.fst.$pop1.$pop2.int $BIN/angsd.Fst.spini.sh $nt $fst_dir $angsd_dir $prefix $pop1 $pop2 no no $denovoref $boot
    for bootv in $( seq 1 $boot ) ; do 
      sbatch --cpus-per-task=${nt} --mem=${memo}G -J ang.Fst.$pop1.$pop2 -o $BIN/temp_out_files/ang.fst.$pop1.$pop2.bt.$bootv $BIN/angsd.Fst.spini.boot.sh $nt $fst_dir $angsd_dir $prefix $pop1 $pop2 no no $denovoref $bootv 
    done
fi ; fi ; fi ; done  
squeue -u $user    
# gather Fst data:
module load system/R-3.5.2 ; Rscript $BIN/gather_fst.R $DATA/sampling.sites $prefix.spini_all $fst_dir $boot

######### individual SFS for HE estimates in angsd #####
prefix="spini"
hedir=$RES/he ; mkdir $hedir
cat $DATA/data.info | tail -n +2 | cut -f 1 | tail -n +2 | while read subdata ; do echo $subdata
    sample_file=indv
    ls $sb_dir/${subdata}_rmdups_noMTnoCP_${prefix}_sorted.bam > $sb_dir/$prefix.$subdata.${sample_file}.bamlist  # create bam list
    cat $sb_dir/$prefix.$subdata.${sample_file}.bamlist
    nt=4 ; memo=10 ; N=1 ; gmin=na ; gmax=na ; mindepthind=1 ; minind=1
    maxdepthind=$(grep "$subdata"  $cov_dir/ind_cov.$prefix.txt | cut -d " " -f3 | sort -n | tail -1) ; echo "maxdepthind = " $maxdepthind
    echo $nt $sb_dir $genome $angsd_dir $prefix $sample_file $gmin $gmax $maxdepthind $mindepthind $minind $subdata
#estimate genotype likelihood for sfs
    sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ang.he.$prefix.$subdata -o $BIN/temp_out_files/ang.he.ind.$prefix.$subdata $BIN/angsd.he.ind.spini.sh $nt $sb_dir $genome $angsd_dir $prefix $sample_file $gmin $gmax $maxdepthind $mindepthind $minind $subdata $sb_dir/$prefix.$subdata.${sample_file}.bamlist
done
squeue -u $user


########### genetic distances with ngstools (reuse beagle file) ############
ngdir=$RES/ngstools ; mkdir $ngdir
sample_file="spini_all"
prefix="spini"
ls -lh $angsd_dir 
memo=4 ; nt=4
beaglefile=$angsd_dir/$prefix.$sample_file.beagle.gz
outfile=$ngdir/$prefix.$sample_file
N=$(cat $DATA/data.info | tail -n +2 | wc -l | cut -d " " -f1)
sbatch --cpus-per-task=${nt} --mem=${memo}G -J ads.$prefix.$sample_file -o $BIN/ads.$prefix.$sample_file $BIN/angsd_dist.sh $beaglefile $outfile $nt $DATA $angsd_dir $prefix $sample_file $N $BIN
squeue -u $user 
ls -lh $ngdir/$prefix.${sample_file}.*


########### scrap ########


awk '{print $6
for loc in xchrom ; do #autosomes
  prefix=rad.ptri
  order_file=$data/pango_sample_list.txt
  sample_file=pango.rad.all_$loc
  ngsadmix_dir=$res/ngsadmix.$loc
  module load system/R-3.4.3
  Rscript --vanilla $bin/ngsadmix_plots.R $popmap2 $prefix.${sample_file} $ngsadmix_dir $order_file
  module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3
  Rscript --vanilla $bin/ngsadmix_geoplots.R $popmap2 $prefix.${sample_file} $ngsadmix_dir $order_file
done

##for lowread autosomes
loc=autosomes
dir=$res/ngsadmix.${loc}.lowreads
covtable=$olddata/RAD_pango_recap.txt
popmap=$dir/highreadsamples.txt
popmap2=$dir/sampledata${loc}.txt ; rm $popmap2; touch $popmap2
#samples=$res/cleanedsamples.txt
awk '{print $1}' $popmap | while read indv; do #$samples
x=$(grep -w $indv $popmap)
y=$(grep -w $indv $covtable)
printf "$x $y\n" | awk '{print $1, $40, $21, $22, $24, $31, $7, $8, $29, $33, $41}' | awk '{$12 = $11 / $5}1' | sed 's/ /\t/g' | sed 's/,//g' >> $popmap2 #$res/ngsadmix.${loc}/sampledata${loc}.txt  #40 for locality, 27 for popn
done  
