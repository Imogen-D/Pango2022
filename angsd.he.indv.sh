#!/bin/sh
#$ -q workq

###before running, make sure batch output name is correct


user=idumville # change to my user name
DATA=/work/$user/pangolins/RAD/data
RES=/work/$user/pangolins/RAD/results
BIN=/work/$user/pangolins/RAD/bin

nt=$1
sb_dir=$2 #this is actually hedir
genome=$3
angsd_dir=$4
prefix=$5
sample_file=$6
gmin=$7
gmax=$8
maxdepthind=$9
mindepthind=${10}
minind=${11}
subdata=${12}
bamlist=${13}

echo $nt $sb_dir $genome $angsd_dir $prefix $sample_file $gmin $gmax $maxdepthind $mindepthind $minind $subdata
# test: # nt=2 ; sb_dir= ; genome=/work/lchikhi/genomes/GCF_000165445.2_Mmur_3.0_genomic.fna ; angsd_dir=/work/lchikhi/RAD/micro_sp3_phylo/results/angsd ; prefix=Msp3 ; subdata=sp3 ; sample_file=phylo_${subdata}
# prefix=spini ; subdata=GB_236_2018 ; sample_file=indv


echo "bamlist ${bamlist} "

module load bioinfo/angsd0.920
module load bioinfo/angsd-0.923

echo "gmin = " $gmin
echo "gmax = " $gmax
echo "minind = " $minind
echo "maxdepthind = " $maxdepthind
echo "mindepthind = " $mindepthind
echo "genome =" $genome
echo "bamfile =" $sb_dir/$prefix.$subdata.${sample_file}.bamlist
echo "nt =" $nt
echo "prefix=" $prefix
echo "subdata =" $subdata
echo "sample file =" $sample_file
N=1

echo "#################################################"
echo "1: run angsd -dosaf"
echo "snpPval not used since we want all sites"
echo "estimate genotype likelihood for unfolded sfs"
angsd -doCounts 1 -setMinDepth $gmin -setMaxDepthInd $maxdepthind -setMinDepthInd $mindepthind -minInd $minind -GL 1 -out $angsd_dir/$prefix.$subdata.${sample_file}.uf -ref $genome -nThreads $nt -doMaf 1 -doSaf 1 -bam $bamlist -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -C 50 -baq 1 -skipTriallelic 1 -only_proper_pairs 1 -anc $genome -doMajorMinor 1 > $angsd_dir/$prefix.$subdata.${sample_file}.out 
echo "1b: moving verbose outfile to new position"
cp $BIN/batch_outputs/ang.he.ind.$prefix.$subdata.$loc $angsd_dir/$prefix.$subdata.${sample_file}.uf.out1 
nline=$(cat $BIN/batch_outputs/ang.he.ind.$prefix.$subdata.$loc | wc -l | cut -d " " -f1)

echo "#################################################"
echo "2: run angsd/realSFS"
echo "do I need maxiter and boostrap ?"
realSFS $angsd_dir/$prefix.$subdata.${sample_file}.uf.saf.idx -P $nt > $angsd_dir/$prefix.$subdata.${sample_file}.uf.sfs
 
echo "#################################################"
echo "3: run angsd -dosaf with -pest .sfs(prior) estmated at the previous steps"
angsd -doCounts 1 -setMinDepth $gmin -setMaxDepthInd $maxdepthind -setMinDepthInd $mindepthind -minInd $minind -GL 1 -out $angsd_dir/$prefix.$subdata.${sample_file}.ufp -ref $genome -nThreads $nt -doMaf 1 -doSaf 1 -doMajorMinor 1 -bam $bamlist -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -C 50 -baq 1 -skipTriallelic 1 -only_proper_pairs 1 -anc $genome -pest $angsd_dir/$prefix.$subdata.${sample_file}.uf.sfs > $angsd_dir/$prefix.$subdata.${sample_file}.ufp.out 
tail -n +$nline $BIN/batch_outputs/ang.he.ind.$prefix.$subdata.$loc > $angsd_dir/$prefix.$subdata.${sample_file}.uf.out2 

echo "#################################################"
echo "4: prepare parameters for ngsStat"
nret_sites2=$( awk '{sum=0; for(i=1; i<=NF; i++) sum += $i; OFMT="%f"; print sum }' $angsd_dir/$prefix.$subdata.${sample_file}.uf.sfs )
echo "there are $nret_sites2 retained sites"

echo "#################################################"
echo "5a: prepare saf file for popGEN/ngsStat"
zcat $angsd_dir/$prefix.$subdata.${sample_file}.ufp.saf.gz > $angsd_dir/$prefix.$subdata.${sample_file}.ufp.saf

echo "5d: run popGEN/ngsStat window analysis"
module load bioinfo/ngsTools-a4d338d
rm $hedir/$prefix.$subdata.${sample_file}.w.pop.stat
ngsStat -npop 1 -postfiles $angsd_dir/$prefix.$subdata.${sample_file}.ufp.saf -nsites $nret_sites2 -iswin 1 -nind $N -block_size 1000 -outfile $sb_dir/$prefix.$subdata.${sample_file}.w.pop.stat
head $sb_dir/$prefix.$subdata.${sample_file}.w.pop.stat
wc -l $sb_dir/$prefix.$subdata.${sample_file}.w.pop.stat
awk '{ total += $4; count++ } END { print total/count }' $sb_dir/$prefix.$subdata.${sample_file}.w.pop.stat > $sb_dir/$prefix.$subdata.${sample_file}.un.he
echo "He:"
awk '{ total += $4; count++ } END { print total/count }' $sb_dir/$prefix.$subdata.${sample_file}.w.pop.stat 

echo "#################################################"
echo "6a: estimate Ho from SFS"
nsites=$( awk '{sum=0; for(i=1; i<=NF; i++) sum += $i; print sum}' $angsd_dir/$prefix.$subdata.${sample_file}.uf.sfs )
nsite=$( cat $angsd_dir/$prefix.$subdata.${sample_file}.uf.sfs | cut -d" " -f 2 )
echo $nsite > $sb_dir/$prefix.$subdata.${sample_file}.un.ho
echo $nsites >> $sb_dir/$prefix.$subdata.${sample_file}.un.ho
echo "scale=6; $nsite/$nsites" | bc >> $sb_dir/$prefix.$subdata.${sample_file}.un.ho
echo "all site: $nsite" 
echo "polymorphic sites: $nsites"
echo "Ho:" 
echo "scale=6; $nsite/$nsites" | bc
