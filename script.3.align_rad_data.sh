#!/bin/sh

# this script does trimmming then aligns all the data to autosomes, X, mito and checks files
 
nt=4 ; memo=10 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i
user=jsalmona
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data
dmx=$data/demultiplexed_data #; mkdir -p $dmx
dmx2=$data/demultiplexed_data_novaseq2 

#
sb_dir=$res/sam-bam ; mkdir -p $sb_dir
cov_dir=$res/coverage ; mkdir -p $cov_dir
gdir=/work/$user/pangolins/ref_genomes
stdir=$res/stacks ; mkdir -p $stdir
#eager_dir=$res/eager ; mkdir -p $eager_dir
genome=$gdir/Jaziri_pseudohap2_scaffolds_HiC
mkdir -p $bin/batch_outputs

cd $bin
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


#############################################
### 0903 REDO novaseq TRIMMING ###
############################################# 
NT=4 ; memo=16 ;
ref=$gdir/Jaziri_pseudohap2_scaffolds_HiC
suffix=_novaseq_1_2021 #2_2021 
reads=pe
prefix="ptri"
save=/save/jsalmona/pangolins/RAD/data
dmxsav=$save/demultiplexed_data #_novaseq2
#dmx2=$data/demultiplexed_data_novaseq2
mkdir -p $res/Trimmomatic_0903_Novaseq1_14 #4-12 #2_1-2 #Trimmomatic_0903_Novaseq2_3-4 #Trimmomatic_0903_Novaseq2_5-9

mkdir -p $bin/batch_outputs_temp

sed -n '334,549 p' $data/RAD_pango_recap.txt > $data/RAD_pango_novaseq2.txt
sed -n '12,333 p' $data/RAD_pango_recap.txt > $data/RAD_pango_novaseq1.txt
#awk -v OFS='\t' '{print $2, $NF}' $data/RAD_pango_novaseq2.txt > $data/RAD_pango_novaseq2_names.txt

for seqs in $(seq 14 14) ; do #$(seq 5 9) 3 4
  echo $ref $prefix $seqs
  indir=$dmxsav.pango${seqs}.bcd.1 #$dmx2.pango${seqs}.bcd.1
  #for files in folder, remove suffix (.1.fq.gz or .2.fq.gz)
  for indv in $(ls $indir/ | sed 's/\..*//' | uniq) ; do
    newname=$(cat $data/RAD_pango_novaseq1.txt | awk -v var="$seqs" '$24 == var' | grep -w $indv | awk '{print $25}')
    #newname=$(grep -w $indv $data/RAD_pango_novaseq2_names.txt | cut -f2)
    F1=$indir/$indv.1.fq.gz
    R2=$indir/$indv.2.fq.gz
    OUT=$res/Trimmomatic_0903_Novaseq1_14 #4-12  #2_1-2 #3-4 #Trimmomatic_0903_Novaseq2_5-9
    echo "#########################################"
    echo "(re)trimming indv ${newname} $F1"
    sbatch --cpus-per-task=${NT} --mem=${memo}G  -J tr.${newname}.pango${seqs} -o $bin/batch_outputs_temp/tr.$prefix.${newname}.$suffix $bin/trimmomatic_low_PE.sh $NT $F1 $R2 $OUT $newname $indir 
    done
  done

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


################################################
### trimming for miseq ###
############################################# 
    
#file organising
sed -n '2,11 p' $data/RAD_pango_recap.txt > $data/RAD_pango_miseq.txt
awk -v OFS='\t' '{print $2, $NF}' $data/RAD_pango_miseq.txt > $data/RAD_pango_miseq_names.txt

ref=$gdir/Jaziri_pseudohap2_scaffolds_HiC
suffix=_miseq_2021
reads=pe
prefix="ptri"

indir=$data/demultiplexed_data.bcd.1
mkdir -p $res/Trimmomatic_0703_Miseq

  #for files in folder, remove suffix (.1.fq.gz or .2.fq.gz)
for indv in $(ls $indir/ | sed 's/\..*//' | uniq) ; do
  newname=$(grep -w $indv $data/RAD_pango_miseq_names.txt | cut -f2)
  F1=$indir/$indv.1.fq.gz
  R2=$indir/$indv.2.fq.gz
  OUT=$res/Trimmomatic_0703_Miseq
  echo "#########################################"
  echo "(re)trimming indv ${newname}"
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J tr.${indv}.pango -o $bin/batch_outputs_temp/tr.$prefix.${indv}.$suffix $bin/trimmomatic_low_PE.sh $NT $F1 $R2 $OUT $newname $indir 
  done
    
  squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"
  
  
for file in * ; do echo $file; done
for file in * ; do head=$(sed 's/_.*//' $file); tail=$(sed 's/\.*//' $file); echo ${head}_1_${tail}; done

  
##############################
#### check all Trimmed there ########
##############################
ls $res/Trimmomatic_* | sed 's/_R.*//' | sort | uniq > temp
awk '{print $NF}' $data/RAD_pango_recap.txt > temp2
cat temp temp2 | sort | uniq -u


##############################
#### index the genome ########
##############################
nt=2 ; memo=10
sbatch --cpus-per-task=${nt} --mem=${memo}G -J index.g -o $bin/batch_outputs/index $bin/index_genome.sh Jaziri_pseudohap2_scaffolds_HiC $gdir 
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"
# check results
ls -lh $gdir
less $bin/batch_outputs/index


###########################################
# aligning trimming miseq to genome#######
###########################################
ref=$gdir/Jaziri_pseudohap2_scaffolds_HiC
NT=4 ; memo=16 ;

  
#file organising
sed -n '2,11 p' $data/RAD_pango_recap.txt > $data/RAD_pango_miseq.txt
awk -v OFS='\t' '{print $2, $NF}' $data/RAD_pango_miseq.txt > $data/RAD_pango_miseq_names.txt

ref=$gdir/Jaziri_pseudohap2_scaffolds_HiC
suffix=_miseq_2021
reads=pe
prefix="ptri"

awk '{print $2}' $data/RAD_pango_miseq_names.txt | while read indv ; do
  indir=$res/Trimmomatic_Miseq
  echo "#########################################"
  echo "processing indv ${indv} with ref $ref in pango miseq"
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J bw.$prefix.${indv}.pango -o $bin/batch_outputs/bw.$prefix.${indv}.$suffix $bin/bwa_miseq2021.sh $NT $indir $indv $sb_dir $ref.fasta $cov_dir $prefix $suffix $reads $memo 
    #$bin/bwa_miseq2021.sh $NT $dmx $indv $sb_dir $ref.fasta $cov_dir $prefix $suffix $reads $memo 
done

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


###########################################
# align novaseq_2021 indv to genome #######
###########################################
ref=$gdir/Jaziri_pseudohap2_scaffolds_HiC
NT=4 ; memo=16 ;
reads=pe
suffix=_novaseq_1_2021 #2_2021 
prefix="ptri"
sb_dir=$res/sam-bam-novaseq1_14 #4-12 #13 #2_3 #1-2 #3-4 #5-9 ; 
mkdir -p $sb_dir
#trim=$res/Trimmomatic_0903_Novaseq2_5-9 #2-3 #1-2

#sed -n '334,549 p' $data/RAD_pango_recap.txt > $data/RAD_pango_novaseq2.txt
#sed -n '12,333 p' $data/RAD_pango_recap.txt > $data/RAD_pango_novaseq1.txt

indir=$res/Trimmomatic_0903_Novaseq1_14 #4-12 #13 #1-2 #2_1-2 #3-4

for indv in $(ls $indir/ | sed 's/_R.*//1' | sort | uniq) ; do
#awk '{print $25}' $data/RAD_pango_novaseq2.txt | while read indv ; do
  echo "#########################################"
  echo "processing indv ${indv} with ref $ref in ${indir}"
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J bw.$prefix.${indv}.pango -o $bin/batch_outputs_temp/bw.$prefix.${indv}.$suffix $bin/bwa_miseq2021.sh $NT $indir $indv $sb_dir $ref.fasta $cov_dir $prefix $suffix $reads $memo 
    #$bin/bwa_miseq2021.sh $NT $dmx $indv $sb_dir $ref.fasta $cov_dir $prefix $suffix $reads $memo 
done

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"


###########################################
# align novaseq_2021 indv to mito genome #######
###########################################
mtref=$gdir/KP306514.1
NT=6 ; memo=30 ;
reads=pe
mtsuffix=_mt_miseq_2021 #_mt_novaseq_1_2021 #2_2021 #_novaseq2021
prefix="ptri"
mtsb_dir=$res/sam-bam-mt-miseq #sam-bam-mt-novaseq1_14 #4-12 #13 #1-2 #2_1-2 #3-4 #5-9 ; 
mkdir -p $mtsb_dir

sed -n '334,549 p' $data/RAD_pango_recap.txt > $data/RAD_pango_novaseq2.txt
sed -n '12,333 p' $data/RAD_pango_recap.txt > $data/RAD_pango_novaseq1.txt

indir=$res/Trimmomatic_Miseq #0903_Novaseq1_14 #4-12 #13 #1-2 #2_1-2 #5-9 #3-4

for indv in $(ls $indir | sed 's/_R.*//1' | sort | uniq) ; do
#awk '{print $25}' $data/RAD_pango_novaseq2.txt | while read indv ; do
  echo "#########################################"
  echo "processing indv ${indv} with ref $mtref in ${indir}"
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J bw.mt.$prefix.${indv}.$mtsuffix -o $bin/batch_outputs_temp/bw.$prefix.${indv}.$mtsuffix $bin/bwa_miseq2021.sh $NT $indir $indv $mtsb_dir $mtref.fasta $cov_dir $prefix $suffix $reads $memo 
    #$bin/bwa_miseq2021.sh $NT $dmx $indv $sb_dir $ref.fasta $cov_dir $prefix $suffix $reads $memo 
done

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

##############################################
#######   check bam files integrity ##########
module load bioinfo/samtools-1.8
#seperately for nova one and two
samtools quickcheck -v $res/sam-bam-novaseq1*/*.rad.ptri.sorted_filtered.bam > $bin/bad_bams1.fofn && echo 'all ok' || echo 'some files failed check, see bad_bams1.fofn'
samtools quickcheck -v $res/sam-bam-novaseq2*/*.rad.ptri.sorted_filtered.bam > $bin/bad_bams2.fofn && echo 'all ok' || echo 'some files failed check, see bad_bams2.fofn'
#seperately for mito alignment
samtools quickcheck -v $res/sam-bam-mt*/*.rad.ptri.sorted_filtered.bam > $bin/bad_bams_mito.fofn && echo 'all ok' || echo 'some files failed check, see bad_bams1.fofn'

grep -i -l "error" $bin/batch_outputs/bw.*
grep -i -l "slurmstepd" $bin/batch_outputs/bw.*
grep -i -l "slurm_script" $bin/batch_outputs/bw.*
grep -i -l "cancel" $bin/batch_outputs/bw.*

##############################################

##############################################
# if some files failed
# 1_check error messages 
grep "slurm_script" $bin/batch_outputs/bw.ptri.*
# 2_rerun failed individuals 
cat $bin/bad_bams1.fofn | cut -d"/" -f8 | cut -d "." -f 1 > $bin/bad_bams1.fofn.ids 
ref=$gdir/Jaziri_pseudohap2_scaffolds_HiC ; NT=6 ; memo=42 ; reads=pe ; suffix=_novaseq_1_2021 ; prefix="ptri"

cat $bin/bad_bams1.fofn.ids | while read indv ; do
  lane=$(ls $res/Trimmomatic_0903_Novaseq1*/${indv}_R1_unpaired*| sed "s|\/${indv}.*||" | sed 's/.*seq1_//')
  indir=$res/Trimmomatic_0903_Novaseq1_${lane}
  echo "#########################################"
  sb_dir=$res/sam-bam-novaseq1_${lane}_2
  mkdir -p $sb_dir
  echo "processing indv ${indv} with ref $ref in $indir to $sb_dir"
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J bw.$prefix.${indv}.pango${seqs} -o $bin/batch_outputs_temp/bw.$prefix.${indv}.$suffix $bin/bwa_miseq2021.sh $NT $indir $indv $sb_dir $ref.fasta $cov_dir $prefix $suffix $reads $memo 
    #$bin/bwa_miseq2021.sh $NT $dmx $indv $sb_dir $ref.fasta $cov_dir $prefix $suffix $reads $memo 
done
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

##now for nova2
cat $bin/bad_bams2.fofn | cut -d"/" -f8 | cut -d "." -f 1 > $bin/bad_bams2.fofn.ids 
cat $bin/bad_bams_recheck.fofn | cut -d"/" -f8 | cut -d "." -f 1 > $bin/bad_bams2.fofn.ids2 
ref=$gdir/Jaziri_pseudohap2_scaffolds_HiC ; NT=10  ; memo=60 ; reads=pe ; suffix=_novaseq_2_2021 ; prefix="ptri" #6 #42

cat $bin/bad_bams2.fofn.ids2 | while read indv ; do
  #indv=Y168_1
  lane=$(ls $res/Trimmomatic_0903_Novaseq2*/${indv}_R1_unpaired*| sed "s|\/${indv}.*||" | sed 's/.*seq2_//')
  indir=$res/Trimmomatic_0903_Novaseq2_${lane}
  echo "#########################################"
  sb_dir=$res/sam-bam-novaseq2_${lane}_Y168_1 #2
  mkdir -p $sb_dir
  echo "processing indv ${indv} with ref $ref in $indir to $sb_dir"
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J bw.$prefix.${indv}.pango${seqs} -o $bin/batch_outputs_temp/bw.$prefix.${indv}.$suffix $bin/bwa_miseq2021.sh $NT $indir $indv $sb_dir $ref.fasta $cov_dir $prefix $suffix $reads $memo 
    #$bin/bwa_miseq2021.sh $NT $dmx $indv $sb_dir $ref.fasta $cov_dir $prefix $suffix $reads $memo 
done
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

##now for mitochondria
cat $bin/bad_bams_mito.fofn | cut -d"/" -f8 | cut -d "." -f 1 > $bin/bad_bams_mito.fofn.ids 
mtref=$gdir/KP306514.1 ; NT=6 ; memo=42 ; reads=pe ; suffix=_mt_novaseq_2021 ; prefix="ptri_mt" #this means name has changed between old and new

cat $bin/bad_bams_mito.fofn.ids | while read indv ; do
  lane=$(ls $res/Trimmomatic_0903_Novaseq*/${indv}_R1_unpaired*| sed "s|\/${indv}.*||" | sed 's/.*seq//')
  indir=$res/Trimmomatic_0903_Novaseq${lane}
  echo "#########################################"
  sb_dir=$res/sam-bam-mt-novaseq${lane}_2
  mkdir -p $sb_dir
  echo "processing indv ${indv} with ref $mtref in $indir to $sb_dir"
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J bw.$prefix.${indv}.pango${lane} -o $bin/batch_outputs_temp/bw.$prefix.${indv}.$suffix $bin/bwa_miseq2021.sh $NT $indir $indv $sb_dir $mtref.fasta $cov_dir $prefix $suffix $reads $memo 
    #$bin/bwa_miseq2021.sh $NT $dmx $indv $sb_dir $mtref.fasta $cov_dir $prefix $suffix $reads $memo 
done
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

##recheck!
module load bioinfo/samtools-1.8

samtools quickcheck -v $res/sam-bam-novaseq*_2/*.rad.ptri.sorted_filtered.bam > $bin/bad_bams_recheck.fofn && echo 'all ok' || echo 'some files failed check, see bad_bams1.fofn'
#DlaA30_3, Fb12_1, T-1866_1 ; Y168_1 done seperately due to still running >24hr

samtools quickcheck -v $res/sam-bam-mt*_2/*.rad.ptri_mt.sorted_filtered.bam > $bin/bad_bams_mito_recheck.fofn && echo 'all ok' || echo 'some files failed check, see bad_bams1.fofn'
#fine

##move new into old files to overwrite old files



###merging all alignments from nova 1 / nova 2 into one file
mkdir -p $res/sam-bam-novaseq1-all
mv $res/sam-bam-novaseq1_*/* $res/sam-bam-novaseq1-all/
mkdir -p $res/sam-bam-novaseq2-all
mv $res/sam-bam-novaseq2_*/* $res/sam-bam-novaseq2-all/

mkdir -p $res/sam-bam-novaseq1-mt
mv -R $res/sam-bam-mt-novaseq1_*/* $res/sam-bam-novaseq1-mt/
mkdir -p $res/sam-bam-novaseq2-mt
mv -R $res/sam-bam-mt-novaseq2_*/* $res/sam-bam-novaseq2-mt/

#removing the empty directories
rmdir $res/sam-bam-novaseq1_*
rmdir $res/sam-bam-novaseq2_*
rmdir $res/sam-bam-mt-novaseq1_*
rmdir $res/sam-bam-mt-novaseq2_*

###################################################################
# Remove X chromosome from all bam files    ####################################
###################################################################
NT=4 ; memo=16
run=miseq #1 #2

prefix="ptri"
#array=( 1-2 3 4-12 14) #for nova 1
#array=( 1-2 3-4 5-9) #for nova 2
indir=$res/sam-bam-miseq-all #novaseq${run}-all
xoutdir=$res/sam-bam-miseq-xchrom ; mkdir -p ${xoutdir} #novaseq${run}-xchrom  
aoutdir=$res/sam-bam-miseq-autosomes ; mkdir -p ${aoutdir} #novaseq${run}-autosomes 
suffix=_miseq_2021 #_novaseq_${run}_2021
reads=pe

gdir=/work/$user/pangolins/ref_genomes

sed '2d' $genome.bed > $gdir/autosomes.bed
 

#for element in  "${array[@]}"; do
#indir=$res/temp

for indv in $(ls ${indir}/ | grep "bai") ; do
  name=$(echo $indv | sed 's/\.rad.ptri.sorted_filtered.bam.bai//')
  infile=${name}.rad.ptri.sorted_filtered.bam
  echo "#########################################"
  echo "getting x of indv ${name} in ${indir}/${infile}"
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J sep.x.$prefix.${name}.pango -o $bin/batch_outputs/x_a_stat.$prefix.${name}.$suffix $bin/sepx_a_stats.sh $NT ${indir} $name ${xoutdir} $ref.fasta $cov_dir $prefix $suffix $reads ${aoutdir} $gdir $memo
  done


squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"




##################################
# harvest a few alignement stats #
##################################


#mkdir -p $bin/alignment_outputs


#mv $bin/batch_outputs_temp $bin/alignment_outputs
indir=$bin/alignment_outputs

echo -e "ID\tqfilt_reads\tmapped_reads\tpercent_mapped\tqmap_reads_all\tqmap_reads_X\tqmap_reads_autosome\tqmap_reads_mt" > $res/miseq_novaseq_2021_alignment_stats # headers

suffix=_novaseq_2_2021 #_novaseq_2_2021 #_miseq_2021  #novaseq_1_2021
prefix="ptri"
awk '{print $25}' $data/RAD_pango_novaseq2.txt | while read indv ; do  #miseq #novaseq2
  echo -e -n "$indv\t" >> $res/miseq_novaseq_2021_alignment_stats # ID indidvidu
  #grep "Input Reads:" $bin/batch_outputs/tri.RI.${indv} | cut -d" " -f 3 | paste -sd+ - | bc | tr '\n' '\t' >> $res/miseq_novaseq_2021_alignment_stats # nb of raw reads
  grep "Total reads:" $indir/bw.$prefix.${indv}.${suffix} | head -n 1 | cut -d ":" -f 2 | tr '\n' '\t' | tr -d ' ' >> $res/miseq_novaseq_2021_alignment_stats # nb of qfiltered reads
  grep "Mapped reads:" $indir/bw.$prefix.${indv}.${suffix} | head -n 1 | cut -d ":" -f 2 | tr '\n' '\t' | tr -d ' ' | tr -d '(' | tr -d ')' >> $res/miseq_novaseq_2021_alignment_stats # nb of aligned reads
  grep "Mapped reads:" $indir/bw.$prefix.${indv}.${suffix} | tail -n +2 | cut -d ":" -f 2 | cut -d "(" -f 1 | tr -d ' ' | tr -d '\n'>> $res/miseq_novaseq_2021_alignment_stats # nb of aligned mapqfilt reads on all, x, noxnomt, and mt
#for x
  grep "Mapped reads:" $indir/x_a_stat.ptri.${indv}.${suffix} | tail -n2  | head -n1 | cut -d ":" -f 2 | cut -d "(" -f 1 | tr -d ' ' | tr -d '\n'>> $res/miseq_novaseq_2021_alignment_stats # nb of aligned mapqfilt reads on all, x, noxnomt, and mt

#for autosome
  grep "Mapped reads:" $indir/x_a_stat.ptri.${indv}.${suffix} | tail -n +2 | cut -d ":" -f 2 | cut -d "(" -f 1 | tr -d ' ' | tr -d '\n'>> $res/miseq_novaseq_2021_alignment_stats # nb of aligned mapqfilt reads on all, x, noxnomt, and mt
#for mt
  grep "Mapped reads:" $indir/bw.$prefix.${indv}._mt${suffix}| tail -n +2 | cut -d ":" -f 2 | cut -d "(" -f 1 | tr -d ' ' | tr -d '\n'>> $res/miseq_novaseq_2021_alignment_stats # nb of aligned mapqfilt reads on all, x, noxnomt, and mt
  echo -e -n "\n" >> $res/miseq_novaseq_2021_alignment_stats
done



#mitos = [_mt_miseq_2021, _mt_novaseq_1_2021, _mt_novaseq_2_2021]


##########
Dla37 is missing from some files (including popmap) ... need to be added back

################################################################
##### Extract individual RAD COVERAGE stats           ##########
################################################################
# grab sample list and number of samples
awk '{print $25}' /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt > $data/pango_sample_list.txt
N=$(cat $data/pango_sample_list.txt | wc -l | cut -d " " -f1)
echo number of individuals = $N

###
NT=1
memo=4
cat $data/pango_sample_list.txt | while read indv ; do
  file=$sb_dir/${indv}.rad.ptri.sorted_filtered.bam
  sbatch --cpus-per-task=${NT} --mem=${memo}G -J star.$indv -o $bin/batch_outputs/cov.$indv $bin/coverage.sh $indv $file $cov_dir
done
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"
###

# harvest coverage stats #
echo -e "ID\tnb_F1R2_sbf1_HF_sites\tnb_F1R2_sbf1_HF_reads\tcov_max\tcov_min\tcov_median\tnb_loci_cov_range\tnb_loci_cov1" > $cov_dir/RAD_cov_stats.txt  # headers
cat $data/pango_sample_list.txt |  while read indv ; do
  #echo $indv
  echo -ne $indv"\t" >> $cov_dir/RAD_cov_stats.txt
  grep $indv.sbf1.f1.bamhits $bin/batch_outputs/cov.$indv | cut -d " " -f1 | tr '\n' '\t' >> $cov_dir/RAD_cov_stats.txt
  sed -ne '/total nb of reads/,$ p' $bin/batch_outputs/cov.$indv | head -n 2 | tail -n 1 | tr '\n' '\t' >> $cov_dir/RAD_cov_stats.txt
  grep quantile $bin/batch_outputs/cov.$indv | cut -d " " -f5 | tr '\n' '\t' | tr -d '"' >> $cov_dir/RAD_cov_stats.txt
  grep median $bin/batch_outputs/cov.$indv | cut -d " " -f4 | tr '\n' '\t' | tr -d '"' >> $cov_dir/RAD_cov_stats.txt
  grep 'loci with cov > 2 x and' $bin/batch_outputs/cov.$indv | cut -d " " -f2 | tr '\n' '\t' | tr -d '"' >> $cov_dir/RAD_cov_stats.txt
  grep 'loci with cov=1' $bin/batch_outputs/cov.$indv | cut -d " " -f2 | tr -d '"' >> $cov_dir/RAD_cov_stats.txt
done
less -S $cov_dir/RAD_cov_stats.txt

# produce coverage plots with all individuals
module load system/R-3.4.3
list=$data/RAD_PG2020_sample_list.txt
Rscript ./cov.plot.all.R $list $cov_dir

##############################################
####MOVE ONTO IDUMVILLE ACCOUNT
##############################################
