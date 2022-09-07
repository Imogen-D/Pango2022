#!/bin/sh

nt=$1 ; echo $1
sb_dir=$2 ; echo $2
popmap=$3 ; echo $3
stdir=$4 ; echo $4
suffix=$5 ; echo $5

module load compiler/gcc-7.2.0
module load bioinfo/stacks-2.5

gstacks -I $sb_dir/ -M $popmap -S $suffix -O $stdir -t $nt

