#!/bin/sh

#this is called by script.3.align_rad_data.sh

module load bioinfo/bwa-0.7.17
module load bioinfo/samtools-1.8
module load bioinfo/GATK-3.8-1-0
module load bioinfo/bamtools-2.5.0

nt=$1 ; echo "nt = $1"
data_dir=$2 ; echo "data_dir = $2"
indv=$3 ; echo "indv = $3"
xoutdir=$4 ; echo "xoutdir = $4"
ref=$5 ; echo "ref= = $5"
cov_dir=$6 ; echo "cov_dir = $6"
prefix=$7 ; echo "prefix = $7"
suffix=$8 ; echo "suffix = $8"
reads=$9 ; echo "reads = $9"
memo=${12} ; echo "memo = ${12}"
memi=$(( memo / nt )) ; echo "memi = $memi"
aoutdir=${10} ; echo "aoutdir = ${10}"
gdir=${11} ; echo "gdir = ${11}"

# nt=4; data_dir=/work/pgaubert/pangolins/RAD/results/Trimmomatic/ ;indv=222492 ; sb_dir=/work/pgaubert/pangolins/RAD/results/sam-bam/ ; ref=/work/pgaubert/pangolins/ref_genomes/Jaziri_pseudohap2_scaffolds_HiC ; cov_dir=/work/pgaubert/pangolins/RAD/results/coverage ; prefix=ptri ; suffix=_miseq2021 ; reads=pe ; memo=8 ; memi=2

echo ""
echo "#########################################################################"
echo "samtools view xchrom" # index the final bam file
time samtools view -b ${data_dir}/${indv}.rad.ptri.sorted_filtered.bam HiC_scaffold_2 > ${xoutdir}/${indv}.sorted_filtered_xchrom.bam

echo ""
echo "#########################################################################"
echo "samtools view autosomes" # index the final bam file
time samtools view -b ${data_dir}/${indv}.rad.ptri.sorted_filtered.bam -M -L $gdir/autosomes.bed > ${aoutdir}/${indv}.sorted_filtered_autosomes.bam
  
#do twice for autosome and xchrom
#out_ext=_sorted_filtered_xchrom.bam ; #define the extensions of our files
out_ext2=sorted_filtered_xchrom.bam 
echo $nt $xoutdir $aoutdir $data_dir $indv $genome $cov_dir $suffix $prefix $sp $reads memi = $memi

echo ""
echo "#########################################################################"
echo "index the final xchrom bam file 2" # index the final bam file
time samtools index -b $xoutdir/${indv}.${out_ext2}

echo "samtools flagstat 2"
samtools flagstat $xoutdir/${indv}.${out_ext2}

echo "bamtools stats 2"  
bamtools stats -in $xoutdir/${indv}.${out_ext2}

echo "plot bamstats 2"
samtools stats $xoutdir/${indv}.${out_ext2} > $xoutdir/${indv}.${out_ext2}.stats
mkdir -p $xoutdir/${indv}.${out_ext2}.stats.folder
/usr/local/bioinfo/src/samtools/samtools-1.8/misc/plot-bamstats $xoutdir/${indv}.${out_ext2}.stats -p $xoutdir/${indv}.${out_ext2}.stats.folder/${indv}.${out_ext2}.stats


#echo ""
#echo "#########################################################################"
#echo "removing temporary data"
#rm $xoutdir/${indv}.${out_ext}*
#rm $xoutdir/${indv}.pe.bam
#rm $xoutdir/${indv}.se.bam

############Now for autosomes
#out_ext_auto=_sorted_filtered_autosomes.bam ; #define the extensions of our files
out_ext2_auto=sorted_filtered_autosomes.bam

echo ""
echo "#########################################################################"
echo "index the final autosome bam file 2" # index the final bam file
time samtools index -b $aoutdir/${indv}.${out_ext2_auto}

echo "samtools flagstat 2 auto"
samtools flagstat $aoutdir/${indv}.${out_ext2_auto}

echo "bamtools stats 2 auto"  
bamtools stats -in $aoutdir/${indv}.${out_ext2_auto}

echo "plot bamstats 2 auto"
samtools stats $aoutdir/${indv}.${out_ext2_auto} > $aoutdir/${indv}.${out_ext2_auto}.stats
mkdir -p $aoutdir/${indv}.${out_ext2_auto}.stats.folder
/usr/local/bioinfo/src/samtools/samtools-1.8/misc/plot-bamstats $aoutdir/${indv}.${out_ext2_auto}.stats -p $aoutdir/${indv}.${out_ext2_auto}.stats.folder/${indv}.${out_ext2_auto}.stats


#echo ""
#echo "#########################################################################"
#echo "removing temporary data"
#rm $aoutdir/${indv}.${out_ext_auto}*
#rm $aoutdir/${indv}.pe.bam
#rm $aoutdir/${indv}.se.bam