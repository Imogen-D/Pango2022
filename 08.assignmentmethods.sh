


nt=1 ; memo=1 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i

# define paths
user=idumville
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data

oldres=/work/jsalmona/pangolins/RAD/results
olddata=/work/jsalmona/pangolins/RAD/data
cd $bin

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

###################################
######## assignPOP #################
#############################
ap_dir=$res/assignment/assignPOP
mkdir -p $res/assignment/assignPOP  


rurlist=$data/WCArurallist.txt

##Okay I'm just running this on autosomes
popdir=$ap_dir; mkdir -p $popdir
rm $ap_dir/WCAruralpopmap.txt ; touch $ap_dir/WCAruralpopmap.txt
popmap=$ap_dir/WCAruralpopmap.txt


awk '{print $1}' $rurlist | while read rurindv; do 
pop=$(grep $rurindv $data/all.popmap)
printf "$pop\n" | sed 's/ /\t/g'  >> $ap_dir/WCAruralpopmap.txt
done

rm $ap_dir/WCAurbanpopmap.txt ; touch $ap_dir/WCAurbanpopmap.txt
popmap=$ap_dir/WCAurbanpopmap.txt
while read indv ; do
printf "$indv Test\n" | sed 's/ /\t/g'  >> $ap_dir/WCAurbanpopmap.txt
done < $data/WCAurban.txt

###################################
######## runnning stacks on these bad bois #################
#############################
## making genepop and structure and hoping like heck doesn't run out of space
ap_dir=$res/assignment/assignPOP
prefix=urban #rural #urban
loc=autosomes
var=assignpop
popmap=$ap_dir/WCA${prefix}popmap.txt
stdir=$res/stacks${loc}
popdir=$ap_dir/${prefix}${var} ; mkdir $popdir


nt=1 ; memo=200


r=0.4 ; R=0.4

sbatch --cpus-per-task=${nt} --mem=${memo}G  -J pop.${var}.${prefix} -o $bin/batch_outputs/gpop.${loc}.${var}${prefix} $bin/run.population.sh $nt $popmap $stdir $popdir $r $R
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

###################################
######## running assignpop #################
#############################
nt=10
memo=10
sbatch --cpus-per-task=${nt} --mem=${memo}G  -J assignpop1 -o $bin/batch_outputs/assignpop1 $bin/run.assignment.sh 

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


###################################
######## making rubias #################
#############################
#running BONE on 100GB for 90%
nt=1
memo=100 #rubias 98 only used 3.8mb oops so 1GB for 98, 8GB for 95, 200GB for 90 failed
  #for BONE tring 200G on 90%

#vcf="/work/idumville/pangolins/RAD/results/locator/concatfull.vcf"
vcf="/work/idumville/pangolins/RAD/results/locator/minmac2wRsalllineagesURpopmap_Doualacomb.vcf"
minmac=minmac2
filt=90 #95 #90 #98 #

sbatch -p unlimitq --cpus-per-task=${nt} --mem=${memo}G  -J BONE${filt}${minmac} -o $bin/batch_outputs/BONE${filt}${minmac} $bin/run.assignment.sh mixref num concatfull col $filt $vcf $minmac

#bone uses mode

##when all done remove 
"/work/idumville/pangolins/RAD/results/locator/concatfull.vcfbgzip.gz"
"/work/idumville/pangolins/RAD/results/locator/minmac2wRsalllineagesURpopmap_Doualacomb.vcfbgzip.gz"
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"



#### running on clusters ## Just minmac2 90 and 95%
minmac=minmac2clusters

#cut -f 5- "/work/idumville/pangolins/RAD/results/assignment/rubiasinput90minmac3.txt" | paste rubiasmetaclusters.txt - | sed 's/ /\t/g' > rubiasinput90${minmac}.txt
#cut -f 5- "/work/idumville/pangolins/RAD/results/assignment/rubiasinput95minmac3.txt" | paste rubiasmetaclusters.txt - | sed 's/ /\t/g' > rubiasinput95${minmac}.txt


#vcf="/work/idumville/pangolins/RAD/results/locator/concatfull.vcf"
vcf="/work/idumville/pangolins/RAD/results/locator/minmac2wRsalllineagesURpopmap_Doualacomb.vcf"

filt=95 #95

nt=1
memo=8

sbatch -p workq --cpus-per-task=${nt} --mem=${memo}G  -J BONE${filt}${minmac} -o $bin/batch_outputs/BONE${filt}${minmac} $bin/run.assignment.sh mixref num concatfull col $filt $vcf $minmac
