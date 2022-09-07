#!/bin/sh
#$ -q workq


vcf=$1
metadata=$2
dir=$3
loc=$4
oldloc=$5

echo load modules
#module unload bioinfo/locator-4760945
#module unload system/Miniconda3

module load system/Miniconda3
module load bioinfo/locator-4760945

#module load bioinfo/tabix-0.2.5
#module load bioinfo/vcftools-0.1.15
#module load bioinfo/bcftools-1.3.1
#bgzip -d -c ${vcf}.gz > ${vcf}

#conda create --name locator
#conda activate locator
#git clone https://github.com/kr-colab/locator.git
#cd locator
#pip install -r req.txt
#conda install -c conda-forge scikit-allel

echo first locator to $dir
locator.py --vcf $vcf --sample_data $metadata --out $dir/firstseed --seed 42111293 #first in $base_loc_dir/seeding.concat.allWCA/seednumbers.txt

echo bootstrapping with 200 boots
locator.py --vcf $vcf --sample_data $metadata --out $dir/bootstrap  --bootstrap --nboots 200 --seed 42111293 #first in $base_loc_dir/seeding.concat.allWCA/seednumbers.txt

#rm ${vcf}
  
exit

echo first locator to $loc
mkdir -p $dir/allsamples/
time locator.py --vcf $vcf --sample_data $metadata --out $dir/${loc}

echo just the first locator is done
exit 

#echo window analysis skipped due to zarr error see script
#zarr.errors.ContainsArrayError: path 'samples' contains an array

echo locator with windows to allsampleswindows
vcf_to_zarr.py --vcf $vcf --zarr $dir/allsnp.zarr
mkdir -p $dir/allsampleswindows/
locator.py --zarr $dir/allsnp.zarr --sample_data $metadata --out $dir/allsampleswindows/ --windows --window_size 250000

#quick look at uncertainty
echo quick uncertainty check with jacknife
mkdir -p $dir/jacknife
locator.py --vcf $vcf --sample_data $metadata --out $dir/jacknife/test --jacknife --nboots 20

#a good look at uncertainty
echo bootstrapping with 5 boots
mkdir -p $dir/bootstrap
locator.py --vcf $vcf --sample_data $metadata --out $dir/bootstrap/test --bootstrap --nboots 5


### plotting ###
echo plotting

module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3

echo plotting my script
Rscript /work/idumville/pangolins/RAD/bin/plottinglocator.R $dir allsamples $metadata $oldloc $loc


Rscript /work/idumville/pangolins/RAD/bin/plot_locator.R --infile $dir/ --sample_data $metadata --out $dir/${loc}_alltest_ --map F #does plot it but can't load map


module unload bioinfo/locator-4760945
module unload system/Miniconda3

exit

echo plot predictions and uncertainties for 9 randomly selected individuals and print the locations with peak kernal density and the geographic centroids across jacknife replicates
Rscript /usr/local/bioinfo/src/locator/locator-4760945/scripts/plot_locator.R --infile $dir/bootstrap/test --sample_data $metadata --out $dir/${loc}_bootstrap_test_
Rscript /usr/local/bioinfo/src/locator/locator-4760945/scripts/plot_locator.R --infile $dir/allsampleswindows --sample_data $metadata --out $dir/${loc}_windows_test_
Rscript /usr/local/bioinfo/src/locator/locator-4760945/scripts/plot_locator.R --infile $dir/jacknife/test --sample_data $metadata --out $dir/${loc}_jackknife_
Rscript /usr/local/bioinfo/src/locator/locator-4760945/scripts/plot_locator.R --infile $dir/ --sample_data $metadata --out $dir/${loc}_alltest_
#infile was /allsampleswindows/

#will be slow
#mkdir $dir/bootstrap
#locator.py --vcf $vcf --sample_data $metadata --out $dir/bootstrap/test --bootstrap --nboots 5

#To unload module, unload in order:
