#!/bin/bash
#SBATCH -p workq

trim_path=/usr/local/bioinfo/src/Trimmomatic/Trimmomatic-0.38/
module load bioinfo/Trimmomatic-0.38

NT=$1
F1=$2
R2=$3
OUT=$4
indv=$5
indir=$6



# NT=4 ; F1=/work/pgaubert/pangolins/RAD/data/demultiplexed_data.pango1.bcd.1/CAM002.1.fq.gz ; R2=/work/pgaubert/pangolins/RAD/data/demultiplexed_data.pango1.bcd.1/CAM002.2.fq.gz  ; indv=CAM002 ; OUT=$res/Trimmomatic/ ; mkdir -p $OUT

echo "run trimmomatic for all $indv runs"
#time java -jar $trim_path/trimmomatic.jar SE -threads $NT $indir/$indv.pb_lastbase.fastq.gz $OUT/$indv.fastq.gz ILLUMINACLIP:$trim_path/adapters/TruSeq3-SE.fa:2:30:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 CROP:98
# NB: for rey-iglesia the minlength needs being very low because around 90% of reads have the index adapter in them :-(

time java -jar $trim_path/trimmomatic.jar PE -threads $NT $F1 $R2 $OUT/"$indv"_R1_paired.fastq.gz $OUT/"$indv"_R1_unpaired.fastq.gz $OUT/"$indv"_R2_paired.fastq.gz $OUT/"$indv"_R2_unpaired.fastq.gz ILLUMINACLIP:$trim_path/adapters/TruSeq3-PE.fa:2:20:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 
echo "compare $indv indir and outdir"
zcat $F1 | wc -l
zcat $OUT/"$indv"_R1_paired.fastq.gz | wc -l
zcat $OUT/"$indv"_R1_unpaired.fastq.gz | wc -l

#ILLUMINACLIP: Remove Illumina adapters provided in the TruSeq3-PE.fa file (provided). Initially Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches. These seeds will be extended and clipped if in the case of paired end reads a score of 30 is reached (about 50 bases), or in the case of single ended reads a score of 10, (about 17 bases).
# LEADING: Remove leading low quality or N bases (below quality 3)
# TRAILING: Remove trailing low quality or N bases (below quality 3)
# SLIDINGWINDOW: Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15
# MINLEN: Drop reads which are less than N bases long after these steps
# CROP: 
