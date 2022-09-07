#!/bin/sh
#SBATCH -p workq

fseq=$1
rseq=$2
out_dir=$3
bcode=$4
ddg=$5
bcdist=$6
e1=$7
e2=$8

mkdir -p $out_dir

#module load compiler/gcc-4.9.1
#module load bioinfo/stacks-2.2
module load compiler/gcc-7.2.0
module load bioinfo/stacks-2.5

if [[ $ddg = "no" ]] ; 
then 

cd $out_dir
echo "process radtag considering 1 barcode only"
process_radtags \
-1 $fseq \
-2 $rseq \
-o $out_dir -b $bcode -r -c -q -i 'gzfastq' --retain_header --disable_rad_check \
-e $e1 \
--filter_illumina #\
#--adapter_1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#--adapter_2 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
#--adapter_mm 3

fi

if [[ $ddg = "yes" ]] ; 
then 

cd $out_dir
echo "process radtag considering 2 barcodes"
# for 2bc
process_radtags \
-1 $fseq \
-2 $rseq \
-o $out_dir -b $bcode -r -c -q -i 'gzfastq' --retain_header --disable_rad_check \
--renz_1 $e1 --renz_2 $e2 \
--inline_inline --filter_illumina \
--barcode_dist_2 2 \
--adapter_1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter_2 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
--adapter_mm 3

fi

exit

#>TruSeq3_IndexedAdapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC[ATTCCTA]ATCTCGTATGCCGTCTTCTGCTTG
#>TruSeq3_UniversalAdapter
#AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

# old setting
--adapter_1 ACACTCTTTCCCTACACGACGCTCTTCCGATCT 
--adapter_2 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT 
--adapter_mm 2

