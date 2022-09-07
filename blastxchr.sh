#!/bin/sh

#This is called by script.3b.chromosomes.sh

g_dir=$1

module load bioinfo/ncbi-blast-2.7.1+

#make the database
makeblastdb -in $res/xblast/Jaziri_pseudohap2_scaffolds_HiC.fasta -dbtype nucl

#blast in outformat 6
blastn -db $g_dir/Jaziri_pseudohap2_scaffolds_HiC.fasta -query $g_dir/boxer_XChrom_CM000039_AAEX04000000.fasta -outfmt 6 -max_target_seqs 10 -out $g_dir/xchr_blast_reference.txt