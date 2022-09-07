#!/bin/sh

#this is called by script.3.align_rad_data.sh

module load bioinfo/bwa-0.7.17
module load bioinfo/samtools-1.8
module load bioinfo/GATK-3.8-1-0
module load bioinfo/bamtools-2.5.0

nt=$1 ; echo "nt = $1"
data_dir=$2 ; echo "data directory = $2"
indv=$3 ; echo "indv = $3"
sb_dir=$4 ; echo "sb_dir = $4"
ref=$5 ; echo "ref= = $5"
cov_dir=$6 ; echo "cov_dir = $6"
prefix=$7 ; echo "prefix = $7"
suffix=$8 ; echo "suffix = $8"
reads=$9 ; echo "reads = $9"
memo=${10} ; echo "memo = ${10}"
memi=$(( memo / nt )) ; echo "memi = $memi"

# nt=4; data_dir=/work/pgaubert/pangolins/RAD/results/Trimmomatic/ ;indv=222492 ; sb_dir=/work/pgaubert/pangolins/RAD/results/sam-bam/ ; ref=/work/pgaubert/pangolins/ref_genomes/Jaziri_pseudohap2_scaffolds_HiC ; cov_dir=/work/pgaubert/pangolins/RAD/results/coverage ; prefix=ptri ; suffix=_miseq2021 ; reads=pe ; memo=8 ; memi=2

out_ext=rad.${prefix}.sorted.bam ; #define the extensions of our files
out_ext2=rad.${prefix}.sorted_filtered.bam 
echo $nt $sb_dir $data_dir $indv $genome $cov_dir $suffix $prefix $sp $reads memi = $memi
if [[ "$reads" == "pe" ]] ; then
  echo "data is paired"
  #time bwa mem $ref $data_dir/${indv}.1.fq.gz $data_dir/${indv}.2.fq.gz -t $nt -R '@RG\tID:MMS\tSM:rad1' | samtools view -bhu -@ $nt | samtools sort -@ $nt -m ${memi}G -o $sb_dir/${indv}.${out_ext} -
  echo ""
  echo "#########################################################################"
  echo "aligning paire-end data"
  time bwa mem $ref $data_dir/${indv}_R1_paired.fastq.gz $data_dir/${indv}_R2_paired.fastq.gz -t $nt -R '@RG\tID:MMS\tSM:rad1' | samtools view -bhu -@ $nt | samtools sort -@ $nt -m ${memi}G -o $sb_dir/${indv}.pe.${out_ext} -
  samtools flagstat $sb_dir/${indv}.pe.${out_ext}
  echo ""
  echo "#########################################################################"
  echo "aligning single-end F1 data"
  time bwa mem $ref $data_dir/${indv}_R1_unpaired.fastq.gz -t $nt -R '@RG\tID:MMS\tSM:rad1' | samtools view -bhu -@ $nt | samtools sort -@ $nt -m ${memi}G -o $sb_dir/${indv}.F1se.${out_ext} -
  samtools flagstat $sb_dir/${indv}.F1se.${out_ext}
  echo ""
  echo "#########################################################################"
  echo "aligning single-end R2 data"
  time bwa mem $ref $data_dir/${indv}_R2_unpaired.fastq.gz -t $nt -R '@RG\tID:MMS\tSM:rad1' | samtools view -bhu -@ $nt | samtools sort -@ $nt -m ${memi}G -o $sb_dir/${indv}.R2se.${out_ext} -
  samtools flagstat $sb_dir/${indv}.R2se.${out_ext}
  echo ""
  echo "#########################################################################"
  echo "Merging PE seF1 and seR2 data"
  samtools merge -f $sb_dir/${indv}.${out_ext} $sb_dir/${indv}.pe.${out_ext} $sb_dir/${indv}.F1se.${out_ext} $sb_dir/${indv}.R2se.${out_ext}
  echo ""
  echo "#########################################################################"
  echo "removing temporary data"
  rm $sb_dir/${indv}.pe.${out_ext}
  rm $sb_dir/${indv}.F1se.${out_ext}
  rm $sb_dir/${indv}.R2se.${out_ext}
# -F 4 suppression des reads non mappés
else
  echo "data is single end"
  #time bwa mem $ref $data_dir/${indv}*.R1.fastq.gz -t $nt -R '@RG\tID:MMS\tSM:reseq' | samtools view -bhu -@ $nt | samtools sort -@ $nt -m ${memi}G -o $sb_dir/${indv}.${out_ext} -
fi

echo ""
echo "#########################################################################"
echo "index the final bam file" # index the final bam file
time samtools index -b $sb_dir/${indv}.${out_ext}

echo "samtools flagstat" 
samtools flagstat $sb_dir/${indv}.${out_ext}

echo "bamtools stats"
bamtools stats -in $sb_dir/${indv}.${out_ext}

# I had to filter singletons and proper pairs on 2 different commands since the AS column slides from 14 to 15 in these 2 read types respectiviely !!!!!!!!!!!!
echo ""
echo "#########################################################################"
echo "filter low quality and secondary alignements etc... on singletons"
time samtools view $sb_dir/${indv}.${out_ext} -h -q 60 -F 261 -@ $nt | \
awk 'substr($0,1,1)=="@" || ($7)=="*"' | \
awk 'substr($0,1,1)=="@" || ($9==0)' | \
awk 'substr($0,1,1)=="@" || int(substr($14,6))>110' | \
samtools view -bhu | samtools sort -@ $nt -m ${memi}G -o $sb_dir/${indv}.se1.bam -
echo "samtools flagstat" 
samtools flagstat $sb_dir/${indv}.se1.bam
# -F261 discards unmapped (4) and scondary alignments (256) and discard paired reads (1) 
# the first awk command keeps all SE reads
# the second awk command keeps all SE reads
# the thrid awk command filters reads based on the alignment score (AS) at column 14

echo ""
echo "#########################################################################"
echo "filter low quality and secondary alignements etc... on singletons"
time samtools view $sb_dir/${indv}.${out_ext} -h -q 60 -f 0x8 -F 260 -@ $nt | \
awk 'substr($0,1,1)=="@" || ($7)=="="' | \
awk 'substr($0,1,1)=="@" || ($9>= 120 && $9<=600) || ($9<=-120 && $9>=-600) || ($9==0)' | \
awk 'substr($0,1,1)=="@" || int(substr($14,6))>110' | \
samtools view -bhu | samtools sort -@ $nt -m ${memi}G -o $sb_dir/${indv}.se.bam -
echo "samtools flagstat" 
samtools flagstat $sb_dir/${indv}.se.bam
# -f 0x8 keep only reads with unmapped pair
# -F260 discards unmapped (4) and scondary alignments (256)
# the first awk command keeps all singletons of PE reads aligned on the same chromosome  
# the second awk command filters reads based on the insert size, while keeping header lines and singletons  
# the thrid awk command filters reads based on the alignment score (AS) at column 15


echo ""
echo "#########################################################################"
echo "filter low quality and secondary alignements etc... on proper pairs"
time samtools view $sb_dir/${indv}.${out_ext} -h -q 60 -f 0x2 -F 260 -@ $nt | \
awk 'substr($0,1,1)=="@" || ($7)=="="' | \
awk 'substr($0,1,1)=="@" || ($9>= 120 && $9<=600) || ($9<=-120 && $9>=-600) || ($9==0)' | \
awk 'substr($0,1,1)=="@" || int(substr($15,6))>110' | \
samtools view -bhu | samtools sort -@ $nt -m ${memi}G -o $sb_dir/${indv}.pe.bam -
echo "samtools flagstat" 
samtools flagstat $sb_dir/${indv}.pe.bam
# -f 0x2 keep only reads with mapped pair
# -F260 discards unmapped (4) and scondary alignments (256)
# the first awk command keeps only PE reads aligned on the same chromosome  
# the second awk command filters reads based on the insert size, while keeping header lines and singletons  
# the thrid awk command filters reads based on the alignment score (AS) at column 15
# below is an old non-funtional version of it
# awk '{split($15, subfield, ":"); if(subfield[3]>120) print $0}' | \
echo ""
echo "#########################################################################"
echo "merge sigletons and proper pairs files"
samtools merge -f $sb_dir/${indv}.${out_ext2} $sb_dir/${indv}.se1.bam $sb_dir/${indv}.se.bam $sb_dir/${indv}.pe.bam

echo ""
echo "#########################################################################"
echo "index the final bam file 2" # index the final bam file
time samtools index -b $sb_dir/${indv}.${out_ext2}

echo "samtools flagstat 2"
samtools flagstat $sb_dir/${indv}.${out_ext2}

echo "bamtools stats 2"  
bamtools stats -in $sb_dir/${indv}.${out_ext2}

echo "plot bamstats 2"
samtools stats $sb_dir/${indv}.${out_ext2} > $sb_dir/${indv}.${out_ext2}.stats
mkdir -p $sb_dir/${indv}.${out_ext2}.stats.folder
/usr/local/bioinfo/src/samtools/samtools-1.8/misc/plot-bamstats $sb_dir/${indv}.${out_ext2}.stats -p $sb_dir/${indv}.${out_ext2}.stats.folder/${indv}.${out_ext2}.stats

##Added because plot bamstats doens't seem to be putting it in the right place?? To check / correct
#mv $sb_dir/${indv}.${out_ext2}-* $sb_dir/${indv}.${out_ext2}.stats.folder/

echo ""
echo "#########################################################################"
echo "removing temporary data"
rm $sb_dir/${indv}.${out_ext}*
rm $sb_dir/${indv}.pe.bam
rm $sb_dir/${indv}.se.bam

