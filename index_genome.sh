#!/bin/sh
#SBATCH -p workq


genome=$1
g_dir=$2

cd $g_dir
module load bioinfo/bwa-0.7.17
module load bioinfo/samtools-1.8
module load bioinfo/picard-2.18.2


# index the genome
bwa index $g_dir/$genome.fasta
samtools faidx $g_dir/$genome.fasta
#module load bioinfo/Java8
rm $g_dir/$genome.dict
java -jar -Xmx4g /usr/local/bioinfo/src/picard-tools/picard-2.18.2/picard.jar CreateSequenceDictionary REFERENCE=$g_dir/$genome.fasta OUTPUT=$g_dir/$genome.dict




