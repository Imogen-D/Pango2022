#!/bin/sh
#$ -q workq

# run as :
# cd ~/work/tutuo_rad/bin ; prefix=tchr1 ; for sample_file in tuto_rad ; do for K in $(seq 1 6) ; do for seed in $(seq 1 10) ; do nt=1 ; qsub -pe parallel_smp $nt -l mem=25G,h_vmem=28G -N adm-${K}-${seed} ~/work/tutuo_rad/bin/ngsadmix.sh $nt $prefix $K $seed $sample_file; done ; done ; done

nt=$1
prefix=$2
echo prefix = $prefix
K=$3
seed=$4
sample_file=$5
gmin=$6
gmax=$7
minind=$8
maxdepthind=$9
bin=${10}
ngsadmix_dir=${11}
angsd_dir=${12}
beagle=${13}

echo -e "\n############################"
echo "load angsd module"
module load bioinfo/angsd0.920

echo run NGSadmix from beagle file for K $K iterations $
#NGSADMIX=$bin/NGSadmix

#echo nt $nt prefix $prefix K $K seed $seed sample_file $sample_file gmin $gmin gmax $gmax minind $minind maxdepthind $maxdepthind bin $bin ngsadmixdir $ngsadmix_dir angsddir $angsd_dir

mkdir -p $ngsadmix_dir/$prefix.${sample_file}_k${K} ; 
NGSadmix -minMaf 0.02 -minInd $minind -likes $beagle -K $K -outfiles $ngsadmix_dir/$prefix.${sample_file}_k${K}/$prefix.${sample_file}_k${K}_s${seed} -P $nt -tol 0.000001 -seed $seed >& $ngsadmix_dir/$prefix.${sample_file}_k${K}/$prefix.${sample_file}_k${K}_s${seed}.err 

rm $ngsadmix_dir/$prefix.${sample_file}_k${K}/$prefix.${sample_file}_k${K}_s${seed}.fopt.gz



