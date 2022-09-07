#!/bin/sh

#this calculates and plots coverage, loci stats

nt=4 ; memo=10 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i

user=idumville
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data

oldres=/work/jsalmona/pangolins/RAD/results


all_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-all
auto_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-autosomes
mt_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-mito
x_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-xchrom


cov_dir=$res/coverage ; mkdir -p $cov_dir
stdir=$res/stacks ; mkdir -p $stdir

gdir=/work/jsalmona/pangolins/ref_genomes
genome=$gdir/Jaziri_pseudohap2_scaffolds_HiC

mkdir -p $bin/batch_outputs

cd $bin
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


########################################
####### COVERAGE #######################
########################################
#get lengths of scaffolds 
module load bioinfo/EMBOSS-6.6.0
infoseq -auto -only -name -length -pgc $genome.fasta > $res/genomelengths.tsv
sed -i -e "1d" $res/genomelengths.tsv
sed -i 's/_/\t/g' $res/genomelengths.tsv

NT=2 ; memo=8

loc=mito #autosomes #autosomes #xchrom #mito = .rad.ptri.sorted_filtered.bam


for indv in $(ls ${mt_sb_dir}/*filtered.bam | sed "s/.rad.ptri.sorted_filtered.bam//" | sed 's/^.*\///' ) ; do #could do this with just autosome too but for now all (#_$loc)
  #indv=12B_1.rad.ptri
  echo $indv from ${mt_sb_dir}/${indv} #_${loc}
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J ${indv}.loci.${loc} -o $bin/batch_outputs/coverage.${indv}.${loc} $bin/coverage.sh $indv $loc $cov_dir $bin $mt_sb_dir
done

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

########################################
####### loci plot #######################
########################################
## unsure where this script has gone eek, I'm assuming it just called the loci.plot.R script

NT=2 ; memo=8
for indv in $(ls ${all_sb_dir}/*filtered.bam |sed "s/\.sorted_filtered.bam//" | sed 's/^.*\///' ) ; do #could do this with just autosome too but for now all (#_$loc)
  #indv=12B_1.rad.ptri
  echo $indv from ${all_sb_dir}/${indv}.sorted_filtered.bam using $bin #_${loc} 
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J ${indv}.loci -o $bin/batch_outputs/loci.${indv} $bin/loci.sh $indv $loc $cov_dir $bin $all_sb_dir
done

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


################################################################
##### Extract individual RAD COVERAGE stats           ##########
################################################################
# grab sample list and number of samples
awk '{print $25}' /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt  > $data/pango_sample_list.txt
N=$(cat $data/pango_sample_list.txt | wc -l | cut -d " " -f1)
echo number of individuals = $N

# harvest coverage stats #
loc=mito #xchrom #autosomes
echo -e "ID\tnb_F1R2_HF_sites\tnb_F1R2_HF_reads\tcov_max\tcov_min\tcov_median\tnb_loci_cov_range\tnb_loci_cov1" > $cov_dir/RAD_${loc}_cov_stats.txt  # headers
cat $data/pango_sample_list.txt |  while read indv ; do
  #echo $indv
  #only for xchrom remove 199_1
  #if [[ "$indv" == 'Y199_1' ]]; then
  #  continue
  #fi
  echo -ne "\n"$indv"\t" >> $cov_dir/RAD_${loc}_cov_stats.txt #no need for "\n" when not mito
  grep $indv.rad.ptri.all $bin/batch_outputs/coverage.$indv.$loc | cut -d " " -f1 | tr '\n' '\t' >> $cov_dir/RAD_${loc}_cov_stats.txt
  sed -ne '/total nb of reads/,$ p' $bin/batch_outputs/coverage.$indv.$loc | head -n 2 | tail -n 1 | tr '\n' '\t' >> $cov_dir/RAD_${loc}_cov_stats.txt
  grep quantile $bin/batch_outputs/coverage.$indv.$loc | cut -d " " -f5 | tr '\n' '\t' | tr -d '"' >> $cov_dir/RAD_${loc}_cov_stats.txt
  grep median $bin/batch_outputs/coverage.$indv.$loc | cut -d " " -f4 | tr '\n' '\t' | tr -d '"' >> $cov_dir/RAD_${loc}_cov_stats.txt
  grep 'loci with cov > 2 x and' $bin/batch_outputs/coverage.$indv.$loc | cut -d " " -f2 | tr '\n' '\t' | tr -d '"' >> $cov_dir/RAD_${loc}_cov_stats.txt
  grep 'loci with cov=1' $bin/batch_outputs/coverage.$indv.$loc | cut -d " " -f2 | tr -d '"' >> $cov_dir/RAD_${loc}_cov_stats.txt
done
less -S $cov_dir/RAD_${loc}_cov_stats.txt


#for mito
sed -i '/^$/d' $cov_dir/RAD_${loc}_cov_stats.txt
#lots with no cov
# produce coverage plots with all individuals


#loc=xchrom
module load system/R-3.4.3
#grep -v "Y199_1" $data/pango_sample_list.txt > $data/${loc}pango_sample_list.txt
list=$data/pango_sample_list.txt #${loc}
#have to use different due to different scales
#Rscript ./xcov.plot.all.R $list $cov_dir $loc &
Rscript ./cov.plot.all.R $list $cov_dir $loc &
