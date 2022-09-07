#!/bin/sh
# This script indexes and blasts the mtDNA and Xchromosome as well as seperating the X from autosomes

nt=4 ; memo=10 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i

user=jsalmona
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data
dmx=$data/demultiplexed_data ; mkdir -p $dmx
#
sb_dir=$res/sam-bam ; mkdir -p $sb_dir
cov_dir=$res/coverage ; mkdir -p $cov_dir
gdir=/work/$user/pangolins/ref_genomes
stdir=$res/stacks ; mkdir -p $stdir
genome=$gdir/Jaziri_pseudohap2_scaffolds_HiC
mkdir -p $bin/batch_outputs

cd $bin
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


##############################
#### index the mtDNA genome ########
##############################
nt=2 ; memo=4
sbatch --cpus-per-task=${nt} --mem=${memo}G -J index.g.mt -o $bin/batch_outputs/index.g.mt $bin/index_genome.sh KP306514.1 $gdir 
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"
# check results
ls -lh $gdir
less $bin/batch_outputs/index.g.mt

##############################
#### index the X chrom ########
##############################
nt=2 ; memo=4
sbatch --cpus-per-task=${nt} --mem=${memo}G -J index.g.xchr -o $bin/batch_outputs/index.g.xchr $bin/index_genome.sh boxer_XChrom_CM000039_AAEX04000000 $gdir 
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"
# check results
ls -lh $gdir
less $bin/batch_outputs/index.g.chro


###################################################################
# create bed files ################################################
###################################################################
prefix="ptri"
bedfile=$genome.bed
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $genome.fasta.fai > $bedfile ; 
head -n 120 $bedfile
wc -l $bedfile

mtref=$gdir/KP306514.1.fasta
mtfile=$gdir/$prefix.mt.bed
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $mtref.fasta.fai > $mtfile ; head $mtfile

Xref=$gdir/boxer_XChrom_CM000039_AAEX04000000
Xfile=$gdir/$prefix.xchro.bed
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $Xref.fasta.fai > $Xfile ; head $Xfile

###################################################################
# Blast mtDNA genome to ref    ####################################
###################################################################
nt=4 ; memo=10
cd /work/jsalmona/pangolins/RAD/results; mkdir -p mitoblast; cd mitoblast
user=jsalmona
gdir=/work/$user/pangolins/ref_genomes

sbatch --cpus-per-task=${nt} --mem=${memo}G -J blastmito -o $bin/batch_outputs/blastmito.g.mt $bin/blastmito.sh $gdir
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

#checking
head mito_blast_reference.txt
#outformat 6 is qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

###################################################################
# Blast X chromosome to ref    ####################################
###################################################################
nt=8 ; memo=20
cd /work/jsalmona/pangolins/RAD/results; mkdir -p xblast; cd xblast
user=jsalmona
gdir=/work/$user/pangolins/ref_genomes

sbatch --cpus-per-task=${nt} --mem=${memo}G -J blastxchr -o $bin/batch_outputs/blastxchr.g.xchr $bin/blastxchr.sh $gdir

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

#checking
head xchr_blast_reference.txt
#outformat 6 is qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore



###################################################################
# Remove X chromosome from all bam files    ####################################
###################################################################
NT=4 ; memo=16
run=1 #2
gdir=/work/$user/pangolins/ref_genomes
prefix="ptri"
array=( 1-2 3 4-12 12 14) #for nova 1
array=( 1-2 3-4 5-9) #for nova 2
indir=$res/sam-bam-novaseq${run}_
xoutdir=$res/bam-xchrom-nova${run}  ; mkdir -p ${xoutdir}
aoutdir=$res/bam-autsomes-nova${run} ; mkdir -p ${aoutdir}
suffix=_novaseq_${run}_2021
reads=pe

sed '2d' $genome.bed > $gdir/autosomes.bed

for element in  "${array[@]}"; do
for indv in $(ls ${indir}${element}/ | grep "bai") ; do
  name=$(echo $indv | sed 's/\.rad.ptri.sorted_filtered.bam.bai//')
  infile=${indv}.rad.ptri.sorted_filtered.bam
  echo "#########################################"
  echo "getting x of indv ${indv} in ${indir}${element}/${infile}"
  samtools view -b ${indir}${element}/${infile} HiC_scaffold_2 > ${xoutdir}/${indv}_xchrom.bam
  samtools view -b ${indir}${element}/${infile} -M -L $gdir/autosomes.bed > ${aoutdir}/${indv}_autosomes.bam
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J bw.mt.$prefix.${indv}.pango -o $bin/batch_outputs_temp/bw.$prefix.${indv}.$suffix $bin/sepx_a_stats.sh $NT ${indir}${element} $indv ${xoutdir} $ref.fasta $cov_dir $prefix $suffix $reads $memo ${aoutdir}
  done
done