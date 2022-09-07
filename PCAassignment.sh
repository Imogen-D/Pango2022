
# script for making Beagle and other information for PCA assignment
# Will be rerunning all lineage PCA


nt=1 ; memo=1 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i

# define paths

user=idumville
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data

cd $bin

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

mkdir -p $res/PCAassignment ; Passign=$res/PCAassignment

#########################################
### Concating beagle files ##############
#########################################

## this is wrong to fix on Monday
xbeagle="/work/idumville/pangolins/RAD/results/rad.ptri.pango.rad.all_xchrom.beagle.gz"
abeagle="/work/idumville/pangolins/RAD/results/rad.ptri.pango.rad.all_autosomes.beagle.gz"

# editing abeagle to remove Y199_1 as not in xbeagle so WILL BE REMOVED FROM ALL ANALYSES
# from grep -n Y199_1 "/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/sampledataautosomes.txt" is line 130 so delete column 130 then concatenate

gunzip -c $abeagle | cut --complement -f130 | awk "{ print NF }" | sort | uniq -c
gzip -c > $Passign/concatallbeagle.gz
gunzip -c $xbeagle | sed  '1d' | gzip -c | awk "{ print NF }" | sort | uniq -c
>> $Passign/concatallbeagle.gz


#######################################################################################################
####### run ngsadmix                                                                         ########## 
#######################################################################################################
prefix=alllineages #parameters used to be prefix=rad.ptri
loc=concat
sample_file=pango.rad.${loc}_$prefix
angsd_dir=$res/angsd ; mkdir -p $angsd_dir

ngsadmix_dir=$Passign ; mkdir -p $ngsadmix_dir

maxdepthind=1000 #placeholder value
  gmax=$(cat $cov_dir/ind_cov.${sample_file}.txt | cut -f5 | paste -sd+ | bc) ; echo "gmax = " $gmax
  
  #gmin is essentially just two times no of samples
  gmin=$(cat $cov_dir/ind_cov.${sample_file}.txt | cut -f4 | paste -sd+ | bc) ; echo "gmin = " $gmin

  minind=1 ; echo "minind = " $minind
  mindepthind=1 #for initial run
 
  for K in $(seq 1 10) ; do for seed in $(seq 1 20) ; do 
    nt=16 ; memo=32 ; if [[ "$K" == "1" ]] ; then nt=8 ; memo=20 ; fi #doubled these for autosomes was 4;8 and 1;2 ish     
        echo $nt $prefix $K $seed $sample_file $gmin $gmax $minind $maxdepthind $bin $ngsadmix_dir $angsd_dir ;  
        sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ngad.$sample_file.${K}.$seed -o $bin/batch_outputs/ngad.$sample_file.${K}.$seed.2 $bin/ngsadmix.sh $nt $prefix $K $seed $sample_file $gmin $gmax $minind $maxdepthind $bin $ngsadmix_dir $angsd_dir;
done ; done 
# monitor your jobs :
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

#######################################################
###################### PCANGSD ########################
#######################################################

mkdir -p $res/PCAassignment ; Passign=$res/PCAassignment
loc=concat
prefix=all
fileprefix=alllineages #ruralhighread_wobio #gab

#if not run ngsadmix, manually make $res/ngsadmix.${loc}.$prefix/${loc}${fileprefix}samples.txt ie for allhighread_wobio

covtable=/work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt

popmap2=$Passign/sampledata${loc}${fileprefix}.txt ; rm $popmap2; touch $popmap2

grep -v Y199_1 $covtable | awk '{print $1, $39, $26}' | sed 's/ /\t/g' | sed 's/,//g' >> $popmap2

beaglefile=$Passign/concatallbeagle.gz

pcangsdir=$Passign

outfile=$pcangsdir/$fileprefix.all
infofile=$popmap2

nt=1 ; memo=4

for npc in all; do # 2 4 6 8  for making distance tables for pcangsd
outfile=$pcangsdir/$fileprefix.$npc
  echo $bin/batch_outputs/pc.$npc.$fileprefix $bin/pcangsd_pango.sh $beaglefile $outfile $nt $npc $pcangsdir $bin $infofile $loc
  sbatch --cpus-per-task=${nt} --mem=${memo}G  -J pc.$npc.$fileprefix -o $bin/batch_outputs/pc.$npc.$fileprefix.${loc} $bin/pcangsd_pango.sh $beaglefile $outfile $nt $npc $pcangsdir $bin $infofile $loc 
done
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

###########################################
##### replotting ##########################
###########################################
###Manually edited sampledataautosomesruralhighread_wobio.txt to include number for NGSADMIX colouring

