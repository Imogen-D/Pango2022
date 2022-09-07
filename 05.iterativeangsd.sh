
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




###########################################################################
########## generate unfolded SFS for pogen He measures ##########
########## at the INDIVIDUAL level for only 10 samples           ##########
###########################################################################
prefix=testing_iterative
loc=autosomes
hedir=$res/${loc}_indv_he${prefix} ; mkdir -p $hedir

#manually made testingsamples.txt with the 3 lowest coverage, 4 above the read threshold, and 3 from the middle from average read number

angsd_dir=$res/angsd${prefix} ; mkdir -p $angsd_dir

for gmin in $(seq 2 2 20) ; do #
  prefix_gmin=testing_iterative_$gmin
  cat $hedir/testingsamples.txt | while read subdata ; do
#cat $res/missinghe.txt |  while read subdata ; do
#subdata=DlaA30_3
    sample_file=indv
    rm $hedir/$prefix_gmin.$subdata.${sample_file}.bamlist
    ls $auto_sb_dir/${subdata}.sorted_filtered_${loc}.bam >> $hedir/$prefix_gmin.$subdata.${sample_file}.bamlist
    #ls $sb_dir/${subdata}_rmdups_noMTnoCP_${prefix}_sorted.bam > $sb_dir/$prefix.$subdata.${sample_file}.bamlist  # create bam list
    cat $hedir/$prefix_gmin.$subdata.${sample_file}.bamlist
    #nt=4 ; memo=10 ; N=1 ; gmin=na ; gmax=na ; mindepthind=1 ; minind=1 ; genome=$denovoref
    nt=1 ; memo=4 ; N=1 ; gmax=na ; mindepthind=1 ; minind=1 ; genome=$gdir/Jaziri_pseudohap2_scaffolds_HiC.fasta #2 4
    maxdepthind=$(cat $cov_dir/ind_cov.pango.rad.all_${loc}.txt | cut -f5 | sort -n | tail -1) ; echo "maxdepthind = " $maxdepthind
    echo $nt $auto_sb_dir $genome $angsd_dir $prefix_gmin $sample_file $gmin $gmax $maxdepthind $mindepthind $minind $subdata
#estimate genotype likelihood for sfs
    sbatch --cpus-per-task=${nt} --mem=${memo}G  -J ang.he.$prefix_gmin.$subdata -o $bin/batch_outputs/ang.he.ind.$prefix_gmin.$subdata.$loc $bin/angsd.he.indv.sh $nt $hedir $genome $angsd_dir $prefix_gmin $sample_file $gmin $gmax $maxdepthind $mindepthind $minind $subdata $hedir/$prefix_gmin.$subdata.${sample_file}.bamlist
done
done
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R" 


rm $hedir/${prefix}_${loc}.txt
sample_file=indv
for gmin in  $(seq 2 2 20) ; do #$(seq 2 2 20)
  prefix_gmin=testing_iterative_$gmin
  cat $hedir/testingsamples.txt |  while read indv ; do
  #indv=Y50_2
    he=$(cat $hedir/${prefix_gmin}.${indv}.${sample_file}.un.he)
    echo $he
    #ho=$(cat $hedir/$prefix.$indv.indv.un.ho | sed '3q;d') $res/ngsadmix.${loc}/sampledata${loc}.txt
  x=$(grep -w $indv $res/ngsadmix.${loc}/sampledata${loc}.txt) #popmap2
  echo $x
  sites=$(grep "polymorphic sites:" $bin/batch_outputs/ang.he.ind.${prefix}_${gmin}.$indv.$loc | awk '{print $3}')
  echo $sites
  if [ -z "$sites" ]
  then
  sites=NA
  fi
  if [ -z "$he" ]
  then
  he=NA
  fi
  printf "${x}\t${he}\t${gmin}\t${sites}\n" | sed 's/ /\t/g' >> $hedir/${prefix}_${loc}.txt  #\t${ho}
done
done

#remake the original he
rm $hedir/all_indv_he_${loc}.txt
oldhedir=$res/${loc}_indv_he
prefix="rad.ptri"
cat $data/pango_sample_list.txt |  while read indv ; do
he=$(cat $oldhedir/$prefix.$indv.indv.un.he)
echo $he
#ho=$(cat $hedir/$prefix.$indv.indv.un.ho | sed '3q;d') $res/ngsadmix.${loc}/sampledata${loc}.txt
x=$(grep -w $indv $res/ngsadmix.${loc}/sampledata${loc}.txt) #popmap2
echo $x
sites=$(grep "polymorphic sites:" $bin/batch_outputs/ang.he.ind.rad.ptri.$indv.$loc | awk '{print $3}')
echo $sites
printf "${x}\t${he}\t0\t${sites}\n" | sed 's/ /\t/g' >> $hedir/all_indv_he_${loc}.txt  #\t${ho}
done

cat $hedir/testingsamples.txt |  while read indv ; do grep $indv $hedir/all_indv_he_${loc}.txt >> $hedir/temp1.txt; done
#cut --complement -f 11 $hedir/temp1.txt  > $hedir/temp2.txt
cat $hedir/temp1.txt $hedir/testing_iterative_${loc}.txt > $hedir/allhe_iterations.txt
rm $hedir/temp*

#where did this sampledata${loc}.txt come from ? unsure what some of the values are

#I plotted this locally


