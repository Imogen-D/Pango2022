#!/bin/sh
#SBATCH -p workq

# called by 01.alignmentstats_coverage.sh

echo "################################"
echo "set paths"
indv=$1 # individual 
loc=$2 # bam file
cov_dir=$3 #cov output directory
bin=$4
sb_dir=$5
file=${sb_dir}/${indv}.rad.ptri.sorted_filtered.bam #.sorted_filtered_${loc}.bam #_${loc}

# there is no need to estimate R2F1 coverage since it is very similar to F1R2 but not easily estimated
# therefore I only estimate coverage for HF sbf1 sites 

echo -e '\ngrab coverage for each sbf1 forward strand site'

#getting correct read
#file=/work/idumville/results/test/23B_1.rad.ptri.sorted_filtered.bam
#top = TGCAG
#samtools view -f 0x40 $file | awk '{print $10}' | cut -c-5 | sort | uniq -c | sort -k1 -nr

module load bioinfo/samtools-1.9

samtools view -f 0x40 $file | awk '$10 ~ /^TGCAG/' | cut -f 3,4 | tr '\t' '-' | sort | uniq -c | sort -k1 -nr | sed -E 's/^ *//; s/ /\t/' > $cov_dir/$indv.$loc.f1.bamhits

#remove - and make tab for loci plot
sed 's/-/\t/' $cov_dir/$indv.$loc.f1.bamhits | sed 's/_/\t/g' > $cov_dir/$indv.$loc.f1.bamhits.scaffoldsep

echo -e '\ntotal nb of sbf1 sites'
wc -l $cov_dir/$indv.$loc.f1.bamhits

echo -e '\ntotal nb of reads for these sbf1 sites'
cat $cov_dir/$indv.$loc.f1.bamhits | cut -f 1 | paste -sd+ | bc 

echo -e '\nproduce individual based coverage plots in R'
module load system/R-3.4.3
Rscript $bin/cov.plot.R $indv.$loc $cov_dir 


echo "################################"
echo "call samtools and R modules"
module load bioinfo/samtools-1.8
module load system/R-3.4.3

echo indv $indv loc $loc cov_dir $cov_dir bin $bin sb_dir $sb_dir file 
$file
echo -e '\nproduce individual based loci plots in R'
module load system/R-3.4.3

Rscript $bin/loci.plot.R $indv.$loc $cov_dir

exit

