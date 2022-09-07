#!/bin/sh

#This is called by script.3b.chromosomes.sh

g_dir=$1

module load bioinfo/ncbi-blast-2.7.1+

#make the database
makeblastdb -in $g_dir/Jaziri_pseudohap2_scaffolds_HiC.fasta -dbtype nucl -out /work/jsalmona/pangolins/RAD/results/mitoblast/Jaziri_ncbidb

#blast in outformat 6

blastn -db $g_dir/Jaziri_pseudohap2_scaffolds_HiC.fasta -query $g_dir/KP306514.1.fasta -outfmt 6 -max_target_seqs 10 -out /work/jsalmona/pangolins/RAD/results/mitoblast/mito_blast_reference.txt

