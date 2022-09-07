#!/bin/sh

bin=$1 ; echo $1
indir=$2 ; echo $2
metadata=$3 ; echo $3
prefix=$4 ; echo $4
knownmetadata=$5 ; echo $5


module load system/Miniconda3
module load bioinfo/locator-4760945

mkdir -p $indir/correctpredlocs

##changing longlat THIS WORKS YASSS
for f in $indir/*_predlocs.txt; do echo "Processing $f file.."
awk '{FS="," ; OFS=","; print $2, $1, $3}' $f > ${f}.opp
sed -i 's/,x,y,sampleID,/x,y,sampleID/' ${f}.opp
mv ${f}.opp $indir/correctpredlocs/
done

awk '{OFS="\t"; print $1, $3, $2}' $metadata > ${metadata}.opp
sed -i 's/sampleID\ty\tx/sampleID\tx\ty/' ${metadata}.opp

awk '{OFS="\t"; print $1, $3, $2}' $knownmetadata > ${knownmetadata}.opp
sed -i 's/sampleID\ty\tx/sampleID\tx\ty/' ${knownmetadata}.opp

rm $indir/${prefix}_areas.txt
#touch $indir/${prefix}_areas.txt
rm $indir/${prefix}error_areas.txt ; touch $indir/${prefix}error_areas.txt

Rscript $bin/plot_locator2.R --infile $indir/correctpredlocs --nsamples 1 --sample_data ${knownmetadata}.opp --out $indir/${prefix}error --map T --error T

sed -i 's/\"//g ; /^x/d ; s/^1 //g'  $indir/${prefix}error_areas.txt
sed -i '1 i\Area Sample Level Prob' $indir/${prefix}error_areas.txt
sed -i 's/[\s]+/\t/g' $indir/${prefix}error_areas.txt

exit 




rm $indir/${prefix}error_areas.txt
rm $indir/${prefix}error_windows.pdf

Rscript $bin/plot_locator2.R --infile $indir/correctpredlocs --nsamples 1 --sample_data ${metadata}.opp --out $indir/${prefix} --map T

sed -i 's/\"//g ; /^x/d ; s/^1 //g'  $indir/${prefix}_areas.txt
sed -i '1 i\Area Sample Level Prob' $indir/${prefix}_areas.txt
sed -i 's/[\s]+/\t/g' $indir/${prefix}_areas.txt

exit 

for f in $indir/*_predlocs.txt; do echo "Processing $f file.."
awk '{FS="," ; OFS=","; print $2, $1, $3}' $f > ${f}.opp
done
sed -i 's/,x,y,sampleID,/x,y,sampleID/' ${f}.opp


Rscript $bin/plot_locator2.R --infile $indir/ --sample_data $indir/metadata/test_sample_data.txt --out $indv_loc_dir/test_ --map F

