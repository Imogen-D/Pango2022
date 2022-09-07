#!/bin/sh
#$ -q workq

nt=$1
genome=$2
res=$3
bin=$4

echo -e "\n############################"
echo "load blast module"
module load bioinfo/ncbi-blast-2.13.0+

echo run all-by-all blast
#makeblastdb -in $genome -dbtype nucl
blastn -db $genome -query $genome -outfmt 6 -out $res/all-vs-all.tsv -num_threads $nt