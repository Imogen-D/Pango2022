#!/bin/sh
#SBATCH -p workq

nt=$1
bamlist=$2
#sb_dir=$2
genome=$3
angsd_dir=$4
prefix=$5
sample_file=$6
gmin=$7
gmax=$8
maxdepthind=$9
mindepthind=${10}
minind=${11}

echo -e "\n############################"
echo "load angsd module"
module load bioinfo/angsd0.920

echo -e "\n############################"
echo "Here come the parameter list"
echo "nt = " $nt 
echo "bamlist = " $bamlist #echo "sb_dir = " $sb_dir 
echo "genome = " $genome 
echo "angsd_dir = " $angsd_dir 
echo "prefix = " $prefix 
echo "sample_file = " $sample_file
echo "gmin = " $gmin
echo "gmax = " $gmax
echo "minind = " $minind
echo "maxdepthind = " $maxdepthind
echo "mindepthind = " $mindepthind

echo -e "\n############################"
echo "Run angsd"
angsd -doCounts 1 -setMinDepth $gmin -setMaxDepth $gmax -setMinDepthInd $mindepthind -minInd $minind -GL 1 -out $angsd_dir/$prefix.${sample_file} -ref $genome -nThreads $nt -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam $bamlist -minQ 20 -minMapQ 30 -minMaf 0.02 -uniqueOnly 1 -remove_bads 1 -C 50 -baq 1 -skipTriallelic 1 -only_proper_pairs 1 > $angsd_dir/$prefix.${sample_file}.out

#sb_dir/$prefix.${sample_file}.bamlist
#for now removing -setMaxDepthInd $maxdepthind 


