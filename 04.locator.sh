
######## All accidentally deleted on 20th June becuase the internet sucks ######

## was present: Seeding locator, inital runs, filtering vcf, plotting ahhhh

nt=10 ; memo=10
nt=1 ; memo=1 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i

# define paths

user=idumville
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data

base_loc_dir=$res/locator

cd $bin

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

##############################
## Extracting replicate data #
##############################
prefix=indv #00 #90

rm touch $base_loc_dir/${prefix}repdist.txt ; touch $base_loc_dir/${prefix}repdist.txt
infile=$base_loc_dir/indvidualcentroids.txt # #$base_loc_dir/seeding.concat.allWCA/seedinglocator90/seeding90error_centroids.txt # #
while read indv; do
  x=$(grep ${indv}_ $infile | awk '{print $1, $2, $3}' | sed 's/ /\t/g' | sed 's/"//g' )
  echo $x >> $base_loc_dir/${prefix}repdist.txt 
done < $base_loc_dir/replicatesamples.txt
##manually edited to remove blank lines and add matrix between replicates >2


##############################
## We're going to try bootstrapppinnnggg #
##############################
#on 10 random WCA urban samples

#first make the vcf (again ugh)
nt=1
memo=10
vcffile="/work/idumville/pangolins/RAD/results/locator/concat.allWCA/concatallWCAfiltered.vcf.gz"
oldmetadata="/work/idumville/pangolins/RAD/results/locator/metadata.txt"
rurlist="/work/idumville/pangolins/RAD/data/WCArurallist.txt"
urblist="/work/idumville/pangolins/RAD/data/WCAurban.txt"

N=10

#awk '{print $1}' $urblist | shuf -n $N > $base_loc_dir/randombootssamps.txt

while read indv; do
  indv_loc_dir=$base_loc_dir/${indv}bootstraps ; mkdir -p $indv_loc_dir
  sbatch --cpus-per-task=${nt} --mem=${memo}G -J ${indv}locvcfs -o $indv_loc_dir/${indv}locvcfs_bs $bin/indvlocvcfs.sh $vcffile $indv_loc_dir $indv $oldmetadata $rurlist
done < $base_loc_dir/randombootssamps.txt


nt=1
memo=15
while read indv; do
   indv_loc_dir=$base_loc_dir/${indv}bootstraps ; mkdir -p $indv_loc_dir;
   sbatch --cpus-per-task=${nt} --mem=${memo}G -J ${indv}locboots -o $indv_loc_dir/${indv}loc_bs $bin/locator.sh $indv_loc_dir/${indv}lociindv.vcf.gz $indv_loc_dir/metadata.txt $indv_loc_dir;
done < $base_loc_dir/randombootssamps.txt

rm $base_loc_dir/bootstraparea.txt ; touch $base_loc_dir/bootstraparea.txt 
while read indv; do
   indv_loc_dir=$base_loc_dir/${indv}bootstraps ; mkdir -p $indv_loc_dir;
   sed -n '2,4p' $indv_loc_dir/*areas.txt | awk '{print $0, "BS"}' | sed 's/ /\t/g' >> $base_loc_dir/bootstraparea.txt
   indv_loc_dir=$base_loc_dir/${indv}.allWCA00_lociindv
   sed -n '2,4p' $indv_loc_dir/*areas.txt | awk '{print $0, "Seed"}' | sed 's/ /\t/g' >> $base_loc_dir/bootstraparea.txt
done < $base_loc_dir/randombootssamps.txt
sed -i '1 i\Area Sample Level Prob Run' $base_loc_dir/bootstraparea.txt 
sed -i 's/ /\t/g' $base_loc_dir/bootstraparea.txt

##############################
## We're going to try the lowread samples #
##############################
#on 10 random WCA urban samples

#first make the vcf (again ugh)
nt=1
memo=10
vcffile="/work/idumville/pangolins/RAD/results/locator/concatfull.vcf.gz"
oldmetadata="/work/idumville/pangolins/RAD/results/locator/metadata.txt"
rurlist="/work/idumville/pangolins/RAD/data/WCArurallist.txt"
urblist="/work/idumville/pangolins/RAD/data/WCAurban.txt"

N=10

while read indv; do
  indv_loc_dir=$base_loc_dir/${indv}.allWCA00_lociindv ; mkdir -p $indv_loc_dir
  sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ${indv}locvcfs -o $indv_loc_dir/${indv}locvcfs $bin/indvlocvcfs.sh $vcffile $indv_loc_dir $indv $oldmetadata $rurlist
done < $data/WCAlowurban.txt

#check number of fields matches metdata 
#gunzip -c $res/locator/Y199_1.allWCA00_lociindv/Y199_1lociindv.vcf.gz | head -n 10000 | grep By3_1 | awk --field-separator=" " "{ print NF }"
#minus 9 from this number


prefix=allWCA
loc=concat
nt=4
mem=50
while read indv; do
  indv_loc_dir=$base_loc_dir/${indv}.allWCA00_lociindv ; mkdir -p $indv_loc_dir
  sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ${indv}loclow -o $indv_loc_dir/concat.allWCA.${indv}00_lociindv.log $bin/iterateseedlocator.sh $indv $prefix $res $loc $indv_loc_dir/${indv}lociindv.vcf.gz $indv_loc_dir/metadata.txt $indv_loc_dir
done < $data/WCAlowurban.txt

#plotting low
nt=1
memo=1
prefix=lowread
base_loc_dir=$res/locator

while read indv; do
  indv_loc_dir=$base_loc_dir/${indv}.allWCA00_lociindv ; mkdir -p $indv_loc_dir
  x=$(grep $indv /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | awk '{print $25, $6, $7}' | sed 's/ /\t/g') #make knownmetadata
  sed '$ d' $indv_loc_dir/metadata.txt > $indv_loc_dir/knownmetadata.txt 
  echo $x >> $indv_loc_dir/knownmetadata.txt 
  #sbatch --cpus-per-task=${nt} --mem=${memo}G -J ${indv}locplotting -o $indv_loc_dir/${indv}plotting $bin/plot_locator_sbatch.sh $bin $indv_loc_dir $indv_loc_dir/metadata.txt ${indv}${prefix} $indv_loc_dir/knownmetadata.txt
  #$bin/plot_locator_sbatch.sh $bin $indv_loc_dir $indv_loc_dir/metadata.txt ${indv}${prefix} $indv_loc_dir/knownmetadata.txt > $indv_loc_dir/${indv}plotting 2>&1
done < $data/WCAlowurban.txt

##############################
## We're going to try all da samples on minmac2 #
##############################
#this is not very iterative but now in sarray
nt=4 ; memo=10 #nt=4 and 10gb for locator only w/o vcf production
sbatch --cpus-per-task=${nt} --mem=${memo}G  $bin/script.sarray.locator.sh $data $bin $base_loc_dir #-o $res/locator/minmac2_iterate

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

##checking which ones are done
"/work/idumville/pangolins/RAD/results/locator/Y226_2_minmac2_lociindv/Y226_2.minmac200_90095729_locindv_predlocs.txt"

base_loc_dir=/work/jsalmona/pangolins/RAD/results/templocator/
while read indv; do
  indv_loc_dir=$base_loc_dir/${indv}* #; rm -d $indv_loc_dir
  #mkdir -p $indv_loc_dir
  x=$(ls $indv_loc_dir/*predlocs.txt | wc -l)
  echo $indv $x >> JStemp.txt
done < $data/WCAurban.txt
base_loc_dir=/work/jsalmona/pangolins/RAD/results/templocator/
base_loc_dir=/work/pgaubert/Imogen/pangolins/RAD/results/locator
/work/idumville/pangolins/RAD/data/WCAurban.txt

########################
##### Running 00 and 90% on minmac2 ###
####################################
vcf="/work/idumville/pangolins/RAD/results/locator/minmac2wRsalllineagesURpopmap_Doualacomb.vcfbgzip.gz"
samples="/work/idumville/pangolins/RAD/data/WCAsamples.txt"
prefix=WCAminmac2

loc_dir=$res/locator

rm $loc_dir/${prefix}knownmetadata.txt
rm $loc_dir/${prefix}metadata.txt

head -n 1 $loc_dir/metadata.txt > $loc_dir/${prefix}knownmetadata.txt 
head -n 1 $loc_dir/metadata.txt > $loc_dir/${prefix}metadata.txt 

while read indv; do
  x=$(grep -w $indv /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | awk '{print $25, $6, $7}' | sed 's/ /\t/g') #make knownmetadata
  echo $x >> $loc_dir/${prefix}knownmetadata.txt 
  y=$(grep $indv $loc_dir/metadata.txt)
  echo $y >> $loc_dir/${prefix}metadata.txt 
done < "/work/idumville/pangolins/RAD/data/WCAsamples.txt"

sed -i 's/ /\t/g' $loc_dir/${prefix}metadata.txt
 
al_dir=$res/locator
vcffile="/work/idumville/pangolins/RAD/results/locator/minmac2wRsalllineagesURpopmap_Doualacomb.vcf.gz"
num=00 #00
minmac=minmac3 #minmac3 #change vcf to suit
nt=4 ; memo=50
sbatch --cpus-per-task=${nt} --mem=${memo}G -J ${num}filtvcfs${minmac} -o $bin/batch_outputs/loc${num}locifilt${minmac} $bin/ALvcfmaking.sh $vcffile $num $minmac $al_dir $bin

##############################
## L27_1 and Y360_1 are CA samples so incl CA training #
##############################

CAsamps=( L27_1 Y360_1 )

#cat $data/CAruralist.txt $data/WCArurallist.txt >> $data/CAWCArurallist.txt
rurlist=$data/CAWCArurallist.txt
prefix=CentAfSamples
vcffile="/work/idumville/pangolins/RAD/results/locator/concatfull.vcf.gz" #on the minmac3 dataset cause most analyses
suff=minmac3

for indv in "${CAsamps[@]}"
do
   indv_loc_dir=$base_loc_dir/${indv}_${prefix}_${suff}_lociindv ; mkdir -p $indv_loc_dir
   sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ${indv}${prefix}${suff}locvcfs -o $indv_loc_dir/${indv}${prefix}locvcfs $bin/indvlocvcfs.sh $vcffile $indv_loc_dir $indv $oldmetadata $rurlist
done

loc=concat
nt=4
memo=50
for indv in "${CAsamps[@]}"
do
  indv_loc_dir=$base_loc_dir/${indv}_${prefix}_${suff}_lociindv ; mkdir -p $indv_loc_dir
  sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ${indv}${prefix}${suff}loc -o $indv_loc_dir/${prefix}.${indv}${suff}_lociindv.log $bin/iterateseedlocator.sh $indv $prefix $res $loc $indv_loc_dir/${indv}lociindv.vcf.gz $indv_loc_dir/metadata.txt $indv_loc_dir
done

suff=minmac2
vcffile="/work/idumville/pangolins/RAD/results/locator/minmac2wRsalllineagesURpopmap_Doualacomb.vcf.gz" #on the minmac2 as well why not

for indv in "${CAsamps[@]}"
do
   indv_loc_dir=$base_loc_dir/${indv}_${prefix}_${suff}_lociindv ; mkdir -p $indv_loc_dir
   sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ${indv}${prefix}${suff}locvcfs -o $indv_loc_dir/${indv}${prefix}locvcfs $bin/indvlocvcfs.sh $vcffile $indv_loc_dir $indv $oldmetadata $rurlist
done


loc=concat
nt=4
memo=50
for indv in "${CAsamps[@]}"
do
  indv_loc_dir=$base_loc_dir/${indv}_${prefix}_${suff}_lociindv ; mkdir -p $indv_loc_dir
  sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ${indv}${prefix}${suff}loc -o $indv_loc_dir/${prefix}.${indv}${suff}_lociindv.log $bin/iterateseedlocator.sh $indv $prefix $res $loc $indv_loc_dir/${indv}lociindv.vcf.gz $indv_loc_dir/metadata.txt $indv_loc_dir
done

#plotting the minmac3
CAsamps=( L27_1 Y360_1 )
prefix=CentAfSamples
suff=minmac3
memo=1 ; nt=1
base_loc_dir=$res/locator
recap=
for indv in "${CAsamps[@]}"
do
  indv_loc_dir=$base_loc_dir/${indv}_${prefix}_${suff}_lociindv ; mkdir -p $indv_loc_dir
  x=$(grep $indv /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | awk '{print $25, $6, $7}' | sed 's/ /\t/g')
  sed '$ d' $indv_loc_dir/metadata.txt > $indv_loc_dir/knownmetadata.txt ; echo $x | sed 's/ /\t/g' >> $indv_loc_dir/knownmetadata.txt
  sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ${indv}${prefix}${suff}plot -o $indv_loc_dir/${prefix}.${indv}${suff}_plotting.log $bin/plot_locator_sbatch.sh $bin $indv_loc_dir $indv_loc_dir/metadata.txt $prefix $indv_loc_dir/knownmetadata.txt
done



#########################
### Zipping #######
#########################
nt=1 memo=1
while read indv; do
	sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ${indv}zip -o $bin/batch_outputs/${indv}zip $bin/zipping.sh $res $indv
done < $data/WCArurallist.txt



#########################
### Extracting PA #######
#########################

"/work/idumville/pangolins/RAD/results/stacksautosomes/populationsautosomes.Y222_3minmac3wRs04/populations.log"
stacksautosomes/populationsautosomes.Y222_3minmac3wRs04:
stacksautosomes/populationsautosomes.Y222_3wRs04:
stacksautosomes/populationsconcat.Y222_3concat00seed:
stacksautosomes/populationsconcat.Y222_3concat90seed:
stacksxchrom/populationsxchrom.Y222_3minmac3wRs04
stacksxchrom/populationsxchrom.Y222_3wRs04:

#######################
### Plotting Bootstraps ##
##############
nt=1
memo=1
prefix=bootstrapping
base_loc_dir=$res/locator

while read indv; do
  indv_loc_dir=$base_loc_dir/${indv}bootstraps ; mkdir -p $indv_loc_dir
  mv $indv_loc_dir/bootstrap_bootFULL_predlocs.txt $indv_loc_dir/FULLpredbootstrap.txt #so not invluded in plot
  mv $indv_loc_dir/firstseed_predlocs.txt $indv_loc_dir/firstseedpres.txt #so not invluded in plot
  rm $indv_loc_dir/knownmetadata.txt 
  rm $indv_loc_dir/correctpredlocs/*
  x=$(grep $indv /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | awk '{print $25, $6, $7}' | sed 's/ /\t/g') #make knownmetadata
  sed '$ d' $indv_loc_dir/metadata.txt > $indv_loc_dir/knownmetadata.txt 
  echo $x >> $indv_loc_dir/knownmetadata.txt 
  sbatch --cpus-per-task=${nt} --mem=${memo}G -J ${indv}locvcfs -o $indv_loc_dir/${indv}plotting $bin/plot_locator_sbatch.sh $bin $indv_loc_dir $indv_loc_dir/metadata.txt ${indv}${prefix} $indv_loc_dir/knownmetadata.txt
done < $base_loc_dir/randombootssamps.txt

#then extract error areas and final distnation for compariosn to original

#dammit I overwrote the bootstrap log file

while read indv; do
indv_loc_dir=$base_loc_dir/${indv}.allWCA00_lociindv #indv_loc_dir=$base_loc_dir/${indv}bootstraps
cp $indv_loc_dir/*error_areas.txt $res/locator/moving/
cp $indv_loc_dir/*centroids.txt $res/locator/moving/
cp $indv_loc_dir/*windows.pdf $res/locator/moving/
done < $data/WCAlowurban.txt #$base_loc_dir/randombootssamps.txt


while read indv; do
  indv_loc_dir=$base_loc_dir/${indv}.allWCA00_lociindv ; mkdir -p $indv_loc_dir
  x=$(grep $indv /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | awk '{print $25, $6, $7}' | sed 's/ /\t/g') #make knownmetadata
  sed '$ d' $indv_loc_dir/metadata.txt > $indv_loc_dir/knownmetadata.txt 
  echo $x >> $indv_loc_dir/knownmetadata.txt 
  sbatch --cpus-per-task=${nt} --mem=${memo}G -J ${indv}locplotting -o $indv_loc_dir/${indv}plotting $bin/plot_locator_sbatch.sh $bin $indv_loc_dir $indv_loc_dir/metadata.txt ${indv}${prefix} $indv_loc_dir/knownmetadata.txt
done < $data/WCAlowurban.txt



#################################
### Rerunning locator on all lineages ##
########################################

al_dir=$res/locatorAL ; mkdir -p $al_dir ; cd $al_dir

#making datasets
#grep Rural /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | awk '{print $25, $6, $7}' | sed 's/ /\t/g' > $al_dir/allrural.txt
#rm $al_dir/urban.txt
#grep -v Rural /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | grep -v Yaounde | grep -v Douala | awk '{print $25, $6, $7}' | sed 's/ /\t/g' >> $al_dir/urban.txt
#grep Yaounde /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | awk '{print $25, $6, $7}' | sed 's/ /\t/g' | shuf -n 10 >> $al_dir/urban.txt
#grep Douala /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | awk '{print $25, $6, $7}' | sed 's/ /\t/g' | shuf -n 10 >> $al_dir/urban.txt
#cat $al_dir/urban.txt $al_dir/allrural.txt > $al_dir/knownmetadata.txt
#awk '{print $1, "NA", "NA"}' $al_dir/urban.txt | sed 's/ /\t/g' > $al_dir/metadata.txt ; cat $al_dir/allrural.txt >> $al_dir/metadata.txt
#awk '{print $1}' $al_dir/urban.txt > $al_dir/urbansamples.txt
#awk '{print $1}' $al_dir/allrural.txt > $al_dir/ruralsamples.txt
#awk '{print $1}' $al_dir/metadata.txt > $al_dir/allsamples.txt 
## indiivudal based to run populations and locator

nt=4 ; memo=20 #nt=4 and 10gb for locator
sbatch --cpus-per-task=${nt} --mem=${memo}G  $bin/ALscript.sarray.locator.sh $al_dir $bin $al_dir
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

minmac3 ignored
Y363_1
Y114_1


## Running for 90 and 00%
#Making vcf


al_dir=$res/locatorAL ; mkdir -p $al_dir ; cd $al_dir
#vcffile="/work/idumville/pangolins/RAD/results/locator/minmac2wRsalllineagesURpopmap_Doualacomb.vcf.gz"
vcffile="/work/idumville/pangolins/RAD/results/locator/concatfull.vcf.gz"
num=90 #90
minmac=minmac3 #minmac2 #change vcf to suit
nt=4 ; memo=20
sbatch --cpus-per-task=${nt} --mem=${memo}G -J AL${num}filtvcfs${minmac} -o $bin/batch_outputs/ALvcf${num}locifilt${minmac} $bin/ALvcfmaking.sh $vcffile $num $minmac $al_dir $bin

##still need to run populations on this

#plotting
nt=1
memo=1
base_loc_dir=$res/locatorAL
while read indv; do
  indv_loc_dir=$base_loc_dir/${indv}_ALminmac3_lociindv ; mkdir -p $indv_loc_dir #minmac2 and minmac3
  x=$(grep -w $indv /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | awk '{print $25, $6, $7}' | sed 's/ /\t/g') #make knownmetadata
  sed '$ d' $indv_loc_dir/metadata.txt > $indv_loc_dir/knownmetadata.txt 
  echo $x >> $indv_loc_dir/knownmetadata.txt 
  sbatch --cpus-per-task=${nt} --mem=${memo}G -J ${indv}locplotting -o $indv_loc_dir/${indv}plotting $bin/plot_locator_sbatch.sh $bin $indv_loc_dir $indv_loc_dir/metadata.txt ${indv}${prefix} $indv_loc_dir/knownmetadata.txt
done < $base_loc_dir/urbansamples.txt


nt=1
memo=1
base_loc_dir=$res/locatorAL

while read indv; do
  indv_loc_dir=$base_loc_dir/${indv}_ALminmac2_lociindv
  head -n 4 $indv_loc_dir/*areas.txt > $indv_loc_dir/${indv}error_areas2.txt
  rm $indv_loc_dir/${indv}error_areas.txt ; mv $indv_loc_dir/${indv}error_areas2.txt $indv_loc_dir/${indv}error_areas.txt
done < $base_loc_dir/urbansamples.txt

##plot loci %s 



#plotting summaryREDO WHEN HAVE B51_1 MINMAC2
#rm $base_loc_dir/indvlocmm2predlocs.txt
#rm $base_loc_dir/indvlocmm3predlocs.txt
while read indv; do
  indv_loc_dir=$base_loc_dir/${indv}_ALminmac3_lociindv 
  head -n 2 $indv_loc_dir/*centroids.txt | tail -n 1 >> $base_loc_dir/indvlocmm3predlocs.txt
done < $base_loc_dir/urbansamples.txt

while read indv; do
  grep -w $indv $base_loc_dir/knownmetadata.txt >> $base_loc_dir/unknownsamples.txt
done < $base_loc_dir/urbansamples.txt

sed -i '1 i\sampleID\tx\ty'  $base_loc_dir/unknownsamples.txt

#############################33
######3 extracting PA #############
#################################
st_dir=$res/stacksautosomes

datasets=( concat00seed concat90seed _minmac3 _minmac2 )


rm $st_dir/allPA.txt ; touch $st_dir/allPA.txt
while read indv; do
  for filt in "${datasets[@]}" ; do
  x=$(grep "private alleles: " $st_dir/populationsconcat.${indv}${filt}/*log | awk '{print $NF}')
  echo $indv $filt $x | sed -e  "s/\s\{1,\}/\t/g" >> $st_dir/allPA.txt
  done
done < $data/WCAurban.txt

#### extracting loci in vcf

st_dir=/save/idumville/pangolins/RAD/results/stacksautosomes/

datasets=( concat00seed concat90seed _minmac3 _minmac2 )

while read indv; do
  for filt in "${datasets[@]}" ; do
  x=$(grep -m 2 "Test:" $st_dir/populationsconcat.${indv}${filt}/*log | tail -n 1 | awk '{print $10, $NF}')
  echo $indv $filt $x | sed  "s|/|\t|g" | sed 's/ /\t/g' | sed 's/;//g'  >> $res/locator/allPAandsites.txt
  done
done < $data/WCAurban.txt

sed -i '1 i\sampleID\tfilt\tall\tvariant\tpolymorphic\tPA' $res/locator/allPAandsites.txt

#############################
###### extracting distances ########
#############################
while read indv; do
  Rscript $bin/locator_distances.R /work/idumville/pangolins/RAD/results/locatorAL/${indv}_ALminmac2_lociindv/correctpredlocs $indv ALmm2
done < $res/locatorAL/urbansamples.txt &


while read indv; do
  Rscript $bin/locator_distances.R /work/idumville/pangolins/RAD/results/locator/${indv}.allWCA00_lociindv/correctpredlocs $indv WCAmm3
done < $data/WCAurban.txt &

## for Bootstrapping
while read indv; do
  Rscript $bin/locator_distances.R /work/idumville/pangolins/RAD/results/locator/${indv}bootstraps/correctpredlocs $indv bootstrap
done < $res/locator/randombootssamps.txt &

sed -i 's/$/bs/' "/work/idumville/pangolins/RAD/results/locator/bootstraplocatorsummary.txt"

while read indv; do
  x=$(grep $indv "/work/idumville/pangolins/RAD/results/locator/WCAmm3locatorsummary.txt")
  echo $x "seed" >> "/work/idumville/pangolins/RAD/results/locator/bootstraplocatorsummary.txt"
done < $res/locator/randombootssamps.txt
#mean 0.25 0.5 0.75 max min

## temp to remove, just making vcfs

minmac=minmac3
#vcffile="/work/idumville/pangolins/RAD/results/locator/minmac2wRsalllineagesURpopmap_Doualacomb.vcf.gz"
vcffile="/work/idumville/pangolins/RAD/results/locator/concatfull.vcf.gz"

oldmetadata="/work/idumville/pangolins/RAD/results/locator/metadata.txt"
rurlist="/work/idumville/pangolins/RAD/results/locatorAL/ruralsamples.txt" # 
urblist="/work/idumville/pangolins/RAD/results/locatorAL/urbansamples.txt" # 

minmac2samps=( LL180926-029-A_1 DlaB89_1 )
minmac3samps=( LL180926-022-A_1 LL180926-029-A_1 T-1409_1 DlaB89_1 )

for indv in  "${minmac3samps[@]}"
do
indv_loc_dir=$base_loc_dir/${indv}_${prefix}_lociindv ; mkdir -p $indv_loc_dir
echo $indv_loc_dir
echo vcf production
sbatch --cpus-per-task=1 --mem=50G -J AL${indv}vcf${minmac} -o $bin/batch_outputs/AL${indv}vcf${minmac}  $bin/indvlocvcfs.sh $vcffile /work/idumville/pangolins/RAD/results/locatorAL/${indv}_AL${minmac}_lociindv/ $indv $oldmetadata $rurlist
done

rm $indv_loc_dir/locilist.txt #big file removal


while read indv
do
x=$(ls /work/idumville/pangolins/RAD/results/locatorAL/${indv}*minmac2*/*predlocs* | wc -l)
echo $indv $x
done < "/work/idumville/pangolins/RAD/results/locatorAL/urbansamples.txt"




### plot_locator2 non iterative

infile="/work/idumville/pangolins/RAD/results/locatorAL/B51_1_ALminmac3_lociindv/correctpredlocs"
sample_data="/work/idumville/pangolins/RAD/results/locatorAL/B51_1_ALminmac3_lociindv/knownmetadata.txt.opp"
out="/work/idumville/pangolins/RAD/results/locatorAL/B51_1_ALminmac3_lociindv/test"
width=5
height=4
samples=NULL
nsamples=9
ncol=3
error=T
legend_position="bottom"
map="t"
longlat="T"
map="T"
haploid=FALSE
centroid_method="kd"

############### 
## extracting all predlocs ##
########

base_loc_dir=$res/locatorAL
rm $base_loc_dir/allALpredlocsmm2.txt
touch $base_loc_dir/allALpredlocsmm2.txt

while read indv; do
 indv_loc_dir=$base_loc_dir/${indv}_ALminmac2_lociindv
for file in $indv_loc_dir/*predlocs.txt
do
  tail -n 1 $file | sed 's/,/\t/g' >> $base_loc_dir/allALpredlocsmm2.txt
  done
done < $base_loc_dir/urbansamples.txt
