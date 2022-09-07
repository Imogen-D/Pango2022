#!/bin/sh

nt=1 ; memo=1 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i

user=jsalmona
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data

#

stdir=$res/stacks ; mkdir -p $stdir
mkdir -p $bin/batch_outputs

cd $bin
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

 
#############################################
### run populations ###
############################################# 
#getting the vcf file from another account
#gunzip -c "/work/idumville/pangolins/RAD/results/locator/minmac2wRsalllineagesURpopmap_Doualacomb.vcf.gz" > "/work/jsalmona/pangolins/RAD/results/iterativestacks/allminmac2.vcf"


###Now in sarray
awk '{print $2}' $rawpopmap | sort | uniq > $data/poplist.txt

# CAUTION run.iterativepop needs to be changed depending on minimum number of samples wanted: N2 or N1
nt=1 ; memo=30
sbatch --cpus-per-task=${nt} --mem=${memo}G -o $bin/batch_outputs/gpop.N2.resampling.minmac2 $bin/run.iterativepop.sh $res $data

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


################################
###### Rerun populations on whole dataset with localities instead of popns #######
################################
#this populations.sh is not very iterative (vcf used stated in script)

#awk '{print $25, $39}' $data/RAD_pango_recap.txt | sed '1d' | sed 's/ /\t/g' > metadata.txt
#sed -i 's/Douala_Central_Market/Douala_All_Markets/g' $res/iterativestacks/metadata.txt
#sed -i 's/Douala_Dakat_Market/Douala_All_Markets/g' $res/iterativestacks/metadata.txt
#nt=1 ; memo=100; sbatch --cpus-per-task=${nt} --mem=${memo}G -o $bin/batch_outputs/populationsconcatfull $bin/populations.sh

################################
###### Extract number PA #######
################################

baseoutdir=$res/iterativestacks
rawpopmap=$baseoutdir/popmap.txt
N=1
rm  $res/iterativestacks/${N}resamplingPAsorig.txt ; touch  $res/iterativestacks/${N}resamplingPAsorig.txt
while read pop; do 
awk -v pop="$pop" '$2==pop' $rawpopmap > $baseoutdir/ittemp${N}.txt
  y=$(wc -l $baseoutdir/ittemp${N}.txt | awk '{print $1}')
  if [[ $y -ge $N ]]; then
    a=$(grep -m 1 $pop "/work/jsalmona/pangolins/RAD/results/iterativestacks/concatfull.p.sumstats_summary.tsv" | awk '{print $2}')
    #echo $pop $a
   echo $pop $a original original ${N} | sed 's/ /\t/g' >> $res/iterativestacks/${N}resamplingPAsorig.txt
     fi
   rm $baseoutdir/ittemp${N}.txt
done < $data/poplist.txt

N=1
rm $res/iterativestacks/${N}resamplingPAs.txt ; touch $res/iterativestacks/${N}resamplingPAs.txt
while read pop; do 
awk -v pop="$pop" '$2==pop' $rawpopmap > $baseoutdir/ittemp${N}.txt
  y=$(wc -l $baseoutdir/ittemp${N}.txt | awk '{print $1}')
  if [[ $y -ge $N ]]; then
   for i in {1..50}; do 
    x=$(grep -m 1 $pop $res/iterativestacks/minmac2concat${i}populations${N}/*summary.tsv | awk '{print $2}') 
    echo $pop $x ${i} iterative ${N} | sed 's/ /\t/g' >> $res/iterativestacks/${N}resamplingPAs.txt; done
    fi
   rm $baseoutdir/ittemp${N}.txt
done < $data/poplist.txt

cat $res/iterativestacks/1resamplingPAs.txt $res/iterativestacks/1resamplingPAsorig.txt $res/iterativestacks/2resamplingPAs.txt $res/iterativestacks/2resamplingPAsorig.txt > $baseoutdir/allresamplingPAs.txt
