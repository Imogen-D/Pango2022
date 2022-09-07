#!/bin/sh
#$ -q workq

#this calculates ploidy, contamination, plotting etc

user=idumville
#user=jsalmona
RES=/work/$user/pangolins/RAD/results
BIN=/work/$user/pangolins/RAD/bin
DATA=/work/$user/pangolins/RAD/data

oldres=/work/jsalmona/pangolins/RAD/results
olddata=/work/jsalmona/pangolins/RAD/data

all_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-all
auto_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-autosomes
mt_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-mito
x_sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-xchrom

###Specifyloc here###
loc=autosomes
stdir=$RES/stacks${loc} ; mkdir $stdir
dnsamp=$DATA/denovo_samples ; mkdir $dnsamp
kmsamp=$DATA/kmer_samples ; mkdir $kmsamp
g_dir=/work/jsalmona/pangolins/ref_genomes
genome=$gdir/Jaziri_pseudohap2_scaffolds_HiC.fasta
sb_dir=/work/jsalmona/pangolins/RAD/results/sam-bam-autosomes
cov_dir=$RES/coverage ; mkdir -p $cov_dir
#"/work/idumville/pangolins/RAD/results/coverage/12B_1.autosomes.f1.bamhits"
data_dir=$DATA/denovo_samples #this isn't used in script so placeholder here

cd $BIN

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"
  
####################################################
##### assess ploidy from bam files #################
nquire_dir=$RES/nquire ; mkdir -p $nquire_dir
ref=$genome;  NT=1 ; memo=2 ;  reads="pe" ;suffix=_paired ; prefix="rad.ptri"
cat $DATA/pango_sample_list.txt| while read indv ; do 
  #echo $NT $data_dir $indv $sb_dir $ref $nquire_dir $prefix $suffix $reads $memo
sbatch --cpus-per-task=${NT} --mem=${memo}G -J nQ.${indv} -o $BIN/batch_outputs/nQ.$prefix.${indv} $BIN/nQuire_ploidy.sh $NT $data_dir $indv $sb_dir $ref $nquire_dir $prefix $suffix $reads $memo
done
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"
##### grab outfiles results #################
echo -e "ID\tDiploid_R2\tTriploid_R2\tTetraploid_R2\tDiploid_loglike\tTriploid_loglike\tTetraploid_loglike" > $RES/ploidy_pangolins_all.txt
cat $DATA/pango_sample_list.txt | while read indv ; do
  echo -e -n "$indv\t" >> $RES/ploidy_pangolins_all.txt
  grep "r^2:" $BIN/batch_outputs/nQ.$prefix.${indv} | cut -d ":" -f 2 | tr '\n' '\t' | tr -d ' ' >> $RES/ploidy_pangolins_all.txt
  grep "loglik:" $BIN/batch_outputs/nQ.$prefix.${indv} | head -n 3 | cut -d ":" -f 2 | tr '\n' '\t' | tr -d ' ' >> $RES/ploidy_pangolins_all.txt
  echo -e -n "\n" >> $RES/ploidy_pangolins_all.txt
done


###################################################
##### check for contaminations among samples + plotting ######
###################################################
loc=autosomes ; 
stdir=$RES/stacks${loc} ; mkdir -p $stdir
ddir=$stdir/populations${loc}.wRs04
vdir=$res/vanquish ; mkdir -p $vdir
##### create individual vcf files #####
cat $data/pango_sample_list.txt | while read id ; do 
  echo $id
  memo=4
  echo $ddir $id
  sbatch --mem=${memo}G -J vc.ct.$id -o $bin/batch_outputs/vc.ct.$id $bin/conta_vcf.sh $ddir $id ; done
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"
ls $ddir/*.gz | wc -l

##### run vanquish ##### (change vanquish_ploidy_plot for tongues)
sample_file=$data/pango_sample_list.txt #"/work/idumville/pangolins/RAD/data/tongue_sample_list.txt"
  memo=2
  x=2
  echo $bin/run_vanquish.sh $vdir ${ddir}/pop. $data $sample_file $bin $x
  sbatch --mem=${memo}G -J vqsh.ct -o $bin/batch_outputs/vqsh.ct $bin/run_vanquish.sh $vdir ${ddir}/pop. $data $sample_file $bin $x
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


###################################################
##### I DIDN'T LOOK AT IGV ALIGNMENTS, STOPPED HERE ######
###################################################


# observe a few alignment files with IGV
module load bioinfo/Java8
# prepare a batch file for IGV with only three individuals to observe a RAD loci
echo "new" >$sb_dir/igv.batch
echo "genome /work/lchikhi/genomes/GCF_000165445.2_Mmur_3.0_genomic.fna" >>$sb_dir/igv.batch
  echo "load ../results/sam-bam/A01_2014_rea_rmdups_mmur3.0_sorted.bam" >>$sb_dir/igv.batch
  echo "load ../results/sam-bam/A02_2014_rea_rmdups_mmur3.0_sorted.bam" >>$sb_dir/igv.batch
  echo "snapshotDirectory /work/lchikhi/RAD/micro_sp3_phylo/bin/" >>$sb_dir/igv.batch
  #echo "goto NC_033660.1:152,000-152,700" >>$sb_dir/igv.batch
  echo "goto NC_033660.1:667000-667600" >>$sb_dir/igv.batch
  echo "sort position" >>$sb_dir/igv.batch
  echo "collapse" >>$sb_dir/igv.batch
  echo "snapshot" >>$sb_dir/igv.batch
done
cat $sb_dir/igv.batch 
# run IGV and observe the RADseq alignment
#java -Xmx4G -jar /usr/local/bioinfo/src/IGV/IGV_2.3.11/igv.jar -g $mt_genome
java -Xmx4G -jar /usr/local/bioinfo/src/IGV/IGV_2.3.11/igv.jar -b $sb_dir/igv.batch


############ look at the igv representation and think about the few simple questions in the slides
# you can also run IGV: java -Xmx4G -jar /usr/local/bioinfo/src/IGV/IGV_2.3.11/igv.jar
# and select the files manually

# count the number of reads after alignment and store it in a file in results



# count nb of cleaned aligned reads store that in a file in results
for prefix in "Msp3" ; do
  maindir=~/work/RAD/micro_sp3_phylo/
  nt=1 ; memo=4
  sbatch --cpus-per-task=${nt} --mem=${memo}G  -J al_count.$prefix -o $BIN/al.count.$prefix $BIN/al_read_counts.sh $maindir $prefix ; done
qstat -u $user 



g_dir=$DATA/genome
genome=$g_dir/chr1.fna
mt_genome=$g_dir/mtdna_genome.fa
rm $RES/clean_read_count_trimmo_aligned ; touch $RES/clean_read_count_trimmo_aligned
cat $DATA/data.info | tail -n +2 | cut -f 2 | while read indv ; do
  echo -n -e $indv '\t' >> $RES/clean_read_count_trimmo_aligned
  echo $((`samtools view $RES/sam-bam/"$indv"_chr1_sorted.bam | wc -l` + `samtools view $RES/sam-bam/"$indv"_mt_sorted.bam | wc -l`)) >> $RES/clean_read_count_trimmo_aligned
done

# count the number of reads aligned after removing duplicates and store it in a file in results
g_dir=$DATA/genome
genome=$g_dir/chr1.fna
mt_genome=$g_dir/mtdna_genome.fa
rm $RES/clean_read_count_trimmo_aligned_rmdups ; touch $RES/clean_read_count_trimmo_aligned_rmdups
cat $DATA/data.info | tail -n +2 | cut -f 2 | while read indv ; do
  echo -n -e $indv '\t' >> $RES/clean_read_count_trimmo_aligned_rmdups
  echo $((`samtools view $RES/sam-bam/"$indv"_rmdups_chr1_sorted.bam | wc -l` + `samtools view $RES/sam-bam/"$indv"_rmdups_mt_sorted.bam | wc -l`)) >> $RES/clean_read_count_trimmo_aligned_rmdups
done

# join the different counts of reads
join $RES/raw_read_count $RES/clean_read_count_trimmo | join - $RES/clean_read_count_trimmo_aligned | join - $RES/clean_read_count_trimmo_aligned_rmdups >$RES/read_counts_tot

# plot the different number of reads per individual
Rscript --vanilla $BIN/barplot_nb_reads.R $RES/read_counts_tot