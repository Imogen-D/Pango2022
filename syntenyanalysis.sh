# this script was ran to make the circos input files, which was then plotted on galaxy.eu

nt=1 ; memo=2 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i

user=jsalmona
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data
gdir=/work/$user/pangolins/ref_genomes

genome=$gdir/Jaziri_pseudohap2_scaffolds_HiC
mkdir -p $bin/batch_outputs

cd $bin
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

#############################################
### changing names and extracting top 57 chromosomes of reference ###
############################################# 

cp $gdir/Jaziri_pseudohap2_scaffolds_HiC.fasta $gdir/Jaziri_pseudohap2_scaffolds_HiC_copy.fasta
awk '/^>/{print ">Pango_" ++i; next}{print}' < $gdir/Jaziri_pseudohap2_scaffolds_HiC_copy.fasta > $gdir/Jaziri_renamed.fasta
 
rm $gdir/Jaziri_pseudohap2_scaffolds_HiC_copy.fasta

module load bioinfo/seqkit-v0.16.0
seqkit head -n 57 $gdir/Jaziri_renamed.fasta > $gdir/Jaziri_renamed_57.fasta

##Dog Genome
doggo=$gdir/dog_genome/GCF_000002285.5_Dog10K_Boxer_Tasha_genomic.fna
awk '/^>/{print ">Dog_" ++i; next}{print}' < $doggo > ${doggo}_renamed.fasta
seqkit head -n 39  ${doggo}_renamed.fasta >  ${doggo}_renamed_39.fasta #top 39 from lengths
mv ${doggo}_renamed_39.fasta $gdir/dog_genome/dog_renamed_39.fasta

#cat file to use for karyotype in circos
mum_dir=/work/jsalmona/pangolins/RAD/results/mummer_dog_pango 
cat $gdir/Jaziri_renamed_57.fasta $gdir/dog_genome/dog_renamed_39.fasta > $mum_dir/concatenatedgenomes.fasta


#############################################
## MUMMER alignment ###
############################################# 
nt=1 ; memo=50

out=$res/mummer_dog_pango/ ; mkdir -p $out

sbatch --cpus-per-task=${nt} --mem=${memo}G  -J mummer.pango.dog -o $bin/batch_outputs/mummer.pango_dog_2 $bin/mummer.sh $nt $out $gdir $gdir/Jaziri_renamed_57

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

module load bioinfo/mummer-4.0.0beta2
show-coords -B -c $res/mummer_dog_pango/top_57_pango_dog.delta > $res/mummer_dog_pango/top_57_pango_dog_full_coords.txt 

#keep perc column, only above 70% (also did 90%), get approximate length
awk '{print $1, $9, $10, $8, $11, $12, $13}' ${out}top_57_pango_dog_full_coords.txt   | awk '$7 >= 70'| awk '{ $8 = $6 - $5 } 1' | sed 's/ /\t/g' > ${out}top_57_pango_dog_mummer_abv70.txt

sort -rnk8 ${out}top_57_pango_dog_mummer_abv70.txt | head -n 5000 | awk '{print $1, $2, $3, $4, $5, $6}' | sed 's/ /\t/' > ${out}top5k_pango_dog_abv70.txt

#making lengths of chromosome file for circos
rm $mum_dir/concat_lengths.txt ; touch $mum_dir/concat_lengths.txt
head -n 39 /work/jsalmona/pangolins/ref_genomes/dog_genome/dog_renamed.fasta.fai | awk '{print $1, $2}' | sed 's/ /\t/g' >> $mum_dir/concat_lengths.txt
head -n 57 /work/jsalmona/pangolins/ref_genomes/Jaziri_pseudohap2_scaffolds_HiC.fasta.fai | awk '{print $1, $2}' | sed 's/ /\t/g' >> $mum_dir/concat_lengths.txt

sed -i 's/HiC_scaffold/Pango/' $mum_dir/concat_lengths.txt

