#!/bin/sh

#this is called by syntenyanalysis.sh

nt=$1 ; echo "nt = $1"
out=$2 ; echo "out directory = $2"
gdir=$3 ; echo "gdir = $3"
genome=$4 ; echo "genome = $4"


module load bioinfo/mummer-4.0.0beta2

nucmer --prefix $out/top_57_pango_dog $gdir/dog_genome/dog_renamed_39.fasta $genome.fasta