#Because IBD on angsd script is a mess
#To run IBD, 1) remake list of samples to exclude 2) run beagle filtering 3) run pcagnsd with all 4) concat X + Auto 5) run IBD with all three datasets


nt=1 ; memo=4 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i #--x11 

# define paths
user=idumville
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data
oldres=/work/jsalmona/pangolins/RAD/results
olddata=/work/jsalmona/pangolins/RAD/data
cd $bin
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

IBD_dir=$res/IBD


#######################################################
#########1. move high read samples ####################
#######################################################
## Make list in ngsadmix_plots.R
## High reads on ngsadmix_plots.R (samples after lowest read number with 95% strcuturing at K=6)
#ALL ibd JUST ON CAMEROON
cp  $res/ngsadmix.autosomes/lowreadswoFGGhaWAfCA.txt $IBD_dir #Note this is also removed Gha lineage (two extra samples vs orginal lowread text file)
cp $res/ngsadmix.autosomes/lowreadswoFGGhaWAfCAwoBioko.txt $IBD_dir/lowreadswoFGGhaWAfCAwoBiokoall.txt
cp $res/ngsadmix.autosomes/lowreadswoFGGhaWAfCAwoBiokowoUrban.txt $IBD_dir/lowreadswoFGGhaWAfCAwoBiokorural.txt



#######################################################
######### 2. Filter Beagle ######################
#######################################################

loc=xchrom #autosomes #xchrom
prefix=rural #all #rural

input=$IBD_dir/lowreadswoFGGhaWAfCAwoBioko${prefix}.txt #listtoremove

ngsadmix_dir=/work/idumville/pangolins/RAD/results/ngsadmix.${loc}.lowreads/ ; mkdir -p $ngsadmix_dir

output=${prefix}highread_wobio #ruralhighread_wobio #list to keep #ruralhighread_wobiogab #ruralhighread

nt=1 ; memo=4
sbatch --cpus-per-task=${nt} --mem=${memo}G  -J cutting_beagle${loc}${output} -o $bin/batch_outputs/beaglefilt${loc}${output} $bin/cuttingbeagle.sh ${loc} $res $ngsadmix_dir $output $input
#this outputs to ngsadmix folder so
 
mv ${ngsadmix_dir}${loc}${output}beagle.gz $IBD_dir/
mv $ngsadmix_dir${loc}${output}samples.txt $IBD_dir/

#for concatnated
prefix=all #all #rural
rm $IBD_dir/concat${prefix}highread_wobiobeagle.gz
cp $IBD_dir/autosomes${prefix}highread_wobiobeagle.gz $IBD_dir/concat${prefix}highread_wobiobeagle.gz

gunzip -c $IBD_dir/xchrom${prefix}highread_wobiobeagle.gz | sed  '1d' | gzip -c >> $IBD_dir/concat${prefix}highread_wobiobeagle.gz

# monitor your jobs :
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


# for concatenated data, use loc = autosomes but change beagle and output name $res/ngsadmix.autosomes.lowreads/concatruralhighread_wobiobeagle.gz

#######################################################
######### 3. PCAngsd ######################
#######################################################
### making new infofile
loc=xchrom #autosomes #xchrom
subloc=xchrom #xchrom #concat
prefix=rural #all #rural
samples=$IBD_dir/${loc}${output}samples.txt
output=${prefix}highread_wobio

covtable=$olddata/RAD_pango_recap.txt
popmap=$data/all.popmap
popmap2=$IBD_dir/sampledata${loc}${output}.txt ; rm $popmap2; touch $popmap2
awk '{print $1}' $samples | while read indv; do
x=$(grep -w $indv $popmap)
y=$(grep -w $indv $covtable)
printf "$x $y\n" |  sed 's/ /\t/g' | sed 's/,//g' >> $popmap2
done   

ngsadmix_dir=$res/ngsadmix.$loc.lowreads ; mkdir -p $ngsadmix_dir

beaglefile=$IBD_dir/${subloc}${output}beagle.gz

pcangsdir=$IBD_dir/pcangsd.${output}.${subloc} ; mkdir -p $pcangsdir

outfile=$pcangsdir/$output.${subloc}
infofile=$popmap2 #res/ngsadmix.${loc}/sampledata${loc}.txt

nt=1 ; memo=4 #1,4 for rural and 2, 10 for auto and concat all

npc=all
sbatch --cpus-per-task=${nt} --mem=${memo}G  -J pc.$output.$subloc -o $bin/batch_outputs/pc.$output.$subloc $bin/pcangsd_pango.sh $beaglefile $outfile $nt $npc $pcangsdir $bin $infofile $loc

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

#######################################################
######### 3b. Cutting genetic distance ######################
#######################################################
#filtering beagle isn't working so attempt to use original pcangsd file  $res/pcangsd.xchrom/rad.ptri.pango.rad.all_xchrom.allcov_euc_dist.rda / + autosome
#will cut out columns + rows accroding to urban (WCA+Gab, no seizure); rural (WCA+Gab); cameroon (only Cameroon country) = all high reads (i.e. >15172), $IBD_dir
/lowreadswoFGGhaWAfCAwoBioko${prefix}.txt rad.ptri.pango.rad.all_autosomes.allcov_euc_dist.rda


loc=xchrom #autosomes #xchrom
subloc=xchrom #xchrom #concat
prefix=rural #all #rural

#in R
DGen_all_auto="rad.ptri.pango.rad.all_autosomes.allcov_euc_dist.rda"
DGen_all_auto=load(file = DGen_all_auto)
DGen_all_auto=cov_euc_dd
DGen_all_auto_names <- sub("(.*?)(_)", "", x= rownames(as.matrix(DGen_all_auto))) #becuase atm concatenated
all_names_to_remove <- readLines("../lowreadswoFGGhaWAfCAwoBiokoall.txt")
DGen_all_auto <- as.matrix(DGen_all_auto)
temp <- as.dist(DGen_all_auto[-c(which(DGen_all_auto_names %in% all_names_to_remove)),-c(which(DGen_all_auto_names %in% all_names_to_remove))])


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

