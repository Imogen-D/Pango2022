#!/bin/sh

nt=4 ; memo=10 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i #--x11 

user=jsalmona
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data
dmx=$data/demultiplexed_data ; # mkdir -p $dmx
dmx2=$data/demultiplexed_data_novaseq2 ; # mkdir -p $dmx2

cd $bin
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

################################################
################################################
#### demultiplex data from first MIseq test ####
NT=2 ; memo=4
plate=MIseq_1_2021
fseq=$data/210413_M03493_0323_000000000-JJRBT/Data/Intensities/BaseCalls/PangoA_S1_L001_R1_001.fastq.gz
rseq=$data/210413_M03493_0323_000000000-JJRBT/Data/Intensities/BaseCalls/PangoA_S1_L001_R2_001.fastq.gz
out_dir=$dmx
#bcodefile=$data/test1.miseq.barcodes ; ddr="yes" ; bcdist=3 ; out_dir=$dmx.bcd.$bcdist
bcodefile=$data/test1.miseq.barcode ; ddr="no" ; bcdist=1 ; out_dir=$dmx.bcd.$bcdist
enzyme1=pstI
enzyme2=mspI
sbatch --cpus-per-task=${NT} --mem=${memo}G  -J pcsrd.$plate -o $bin/pcsrd.$bcdist.$plate $bin/process_radtags.pg.sh $fseq $rseq $out_dir $bcodefile $ddr $bcdist $enzyme1 $enzyme2
#$bin/process_radtags.pg.sh $fseq $rseq $out_dir $bcodefile $enzyme1 $enzyme2
squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

#extract relevant information from outfile
head -n 200 $out_dir/process_radtags.BaseCalls.log > process_radtags.BaseCalls.log.bcd.$bcdist.head
################################################

#####################################################
#####################################################
#### demultiplex data from first NOVASEQ plate 1 ####
NT=2 ; memo=8
plate=Novaseq_1_2021
for seqs in $(seq 1 14) ; do # $(seq 1 14) ; do
#seqs=1
  echo "treating rpi pango${seqs}" 
  fseq=$data/RUN1_fasteris/www.fasteris.com/private/ASWZ/ASWZ-20210910/data/210924_A00902_B_L1-2_ASWZ-8-${seqs}_R1.fastq.gz
  rseq=$data/RUN1_fasteris/www.fasteris.com/private/ASWZ/ASWZ-20210910/data/210924_A00902_B_L1-2_ASWZ-8-${seqs}_R2.fastq.gz
  awk -v var="$seqs" '$4 == var' $data/RUN1_novaseq_2021.data.info | awk '{print $7 "\t" $2}' > $data/Novaseq_1_2021.pango${seqs}.barcode 
  bcodefile=$data/Novaseq_1_2021.pango${seqs}.barcode  ; ddr="no" ; bcdist=1 ; out_dir=$dmx.pango${seqs}.bcd.$bcdist # this considers only the first barcode and a barcode distance of 1bp
  enzyme1=pstI ; enzyme2=mspI
  cd $out_dir
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J pcsrd.$plate.$seqs -o $bin/pcsrd.$bcdist.$plate.$seqs $bin/process_radtags.pg.sh $fseq $rseq $out_dir $bcodefile $ddr $bcdist $enzyme1 $enzyme2
  #$bin/process_radtags.pg.sh $fseq $rseq $out_dir $bcodefile $ddr $bcdist $enzyme1 $enzyme2
  cd $bin
done
squeue -u $user -o "%.15i %.35j %.8u %.2t %.10M %.6D %.6C %.6m %R"

#check quickly outfile
plate=Novaseq_1_2021 ; bcdist=1 ; 
for seqs in $(seq 1 14) ; do
  echo "##############################################################"
  echo "$bin/pcsrd.$bcdist.$plate.$seqs"
  tail -n 7 $bin/pcsrd.$bcdist.$plate.$seqs 
  echo "##############################################################"
done
plate=Novaseq_1_2021 ; bcdist=1 ; 


####################################################
# extract relevant information from outfile
# title
seqs=1 ; n1=$( grep -n "NoRadTag" $dmx.pango${seqs}.bcd.$bcdist/process_radtags.data.log | cut -d ":" -f 1 )
head -n $n1 $dmx.pango${seqs}.bcd.$bcdist/process_radtags.data.log | tail -n 1 > $data/demultiplex_summary.all.txt
# miseq run
  echo "##############################################################"
  echo "processing pango${seqs}"
  n1=$( grep -n "NoRadTag" $dmx.bcd.$bcdist/process_radtags.BaseCalls.log | cut -d ":" -f 1 )
  n2=$( grep -n "recorded" $dmx.bcd.$bcdist/process_radtags.BaseCalls.log | cut -d ":" -f 1 )
  head -n $(( $n2 - 2 )) $dmx.bcd.$bcdist/process_radtags.BaseCalls.log | tail -n +$(( $n1 +1)) >> $data/demultiplex_summary.all.txt
# novaseq run 1
for seqs in $(seq 1 14) ; do
  echo "##############################################################"
  echo "processing pango${seqs}"
  n1=$( grep -n "NoRadTag" $dmx.pango${seqs}.bcd.$bcdist/process_radtags.data.log | cut -d ":" -f 1 )
  n2=$( grep -n "recorded" $dmx.pango${seqs}.bcd.$bcdist/process_radtags.data.log | cut -d ":" -f 1 )
  head -n $(( $n2 - 2 )) $dmx.pango${seqs}.bcd.$bcdist/process_radtags.data.log | tail -n +$(( $n1 +1)) >> $data/demultiplex_summary.all.txt
done
head $data/demultiplex_summary.all.txt
wc -l $data/demultiplex_summary.all.txt
####################################################


 #novaseq batch 2; doing 1-5 SEPERATELY as in different folder
  #RUN2_novaseq_2021 data from original barcode file
####www.fasteris.com/private/ASWZ/ASWZ-20211105/data/ now deleted to save space on work folder####
 
######################################################
## demultiplex data from second NOVASEQ plate 2 6-9 ##
######################################################
NT=2 ; memo=8
plate=Novaseq_2_2021
for seqs in $(seq 6 9) ; do 
#seqs=1
  echo "treating rpi pango${seqs} second novaseq" 
  fseq=$data/www.fasteris.com/private/ASWZ/ASWZ-20211105/data/211206_A00902_B_L1-2_ASWZ-9-${seqs}_R1.fastq.gz
  rseq=$data/www.fasteris.com/private/ASWZ/ASWZ-20211105/data/211206_A00902_B_L1-2_ASWZ-9-${seqs}_R2.fastq.gz
  awk -v var="$seqs" '$6 == var' $data/RUN2_novaseq_2021.data.info | awk '{print $3 "\t" $1}' > $data/Novaseq_2_2021.pango${seqs}.barcode 
  bcodefile=$data/Novaseq_2_2021.pango${seqs}.barcode  ; ddr="no" ; bcdist=1 ; out_dir=$dmx2.pango${seqs}.bcd.$bcdist # this considers only the first barcode and a barcode distance of 1bp
  enzyme1=pstI ; enzyme2=mspI
  cd $out_dir
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J pcsrd.$plate.$seqs -o $bin/pcsrd.$bcdist.$plate.$seqs $bin/process_radtags.pg.sh $fseq $rseq $out_dir $bcodefile $ddr $bcdist $enzyme1 $enzyme2
  #$bin/process_radtags.pg.sh $fseq $rseq $out_dir $bcodefile $ddr $bcdist $enzyme1 $enzyme2
  cd $bin
done
squeue -u $user -o "%.15i %.35j %.8u %.2t %.10M %.6D %.6C %.6m %R"


#####################################################
#### demultiplex data from second NOVASEQ plate 2 1-5 ####
#####################################################
##Reran novaseq 1 due to Invalid filename on line 3: 'T-804(EP-4)(Pg-5)' -> manually renamed to T-804 in RUN2_novaseq_2021
NT=2 ; memo=8
plate=Novaseq_2_2021
savepath=/save/jsalmona/pangolins/RAD/data/www.fasteris.com/private/ASWZ/ASWZ-20211105/data
for seqs in $(seq 1 5) ; do  #for seqs in $(seq 1) ; do 
#seqs=1
  echo "treating rpi pango${seqs} second novaseq" 
  fseq=$savepath/211206_A00902_B_L1-2_ASWZ-9-${seqs}_R1.fastq.gz
  rseq=$savepath/211206_A00902_B_L1-2_ASWZ-9-${seqs}_R2.fastq.gz
  awk -v var="$seqs" '$6 == var' $data/RUN2_novaseq_2021.data.info | awk '{print $3 "\t" $1}' > $data/Novaseq_2_2021.pango${seqs}.barcode 
  bcodefile=$data/Novaseq_2_2021.pango${seqs}.barcode  ; ddr="no" ; bcdist=1 ; out_dir=$dmx2.pango${seqs}.bcd.$bcdist # this considers only the first barcode and a barcode distance of 1bp
  enzyme1=pstI ; enzyme2=mspI
  cd $out_dir
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J pcsrd.$plate.$seqs -o $bin/pcsrd.$bcdist.$plate.$seqs $bin/process_radtags.pg.sh $fseq $rseq $out_dir $bcodefile $ddr $bcdist $enzyme1 $enzyme2
  #$bin/process_radtags.pg.sh $fseq $rseq $out_dir $bcodefile $ddr $bcdist $enzyme1 $enzyme2
  cd $bin
done
squeue -u $user -o "%.15i %.35j %.8u %.2t %.10M %.6D %.6C %.6m %R"


#check quickly outfile
plate=Novaseq_2_2021 ; bcdist=1 ; 
for seqs in $(seq 1 9) ; do
  echo "##############################################################"
  echo "$bin/pcsrd.$bcdist.$plate.$seqs"
  tail -n 7 $bin/pcsrd.$bcdist.$plate.$seqs 
  echo "##############################################################"
done


####################################################
# extract relevant information from outfile // appending to summary txt
# novaseq run 2
for seqs in $(seq 1 9) ; do
  echo "##############################################################"
  echo "processing pango${seqs}"
  n1=$( grep -n "NoRadTag" $dmx2.pango${seqs}.bcd.$bcdist/process_radtags.data.log | cut -d ":" -f 1 )
  n2=$( grep -n "recorded" $dmx2.pango${seqs}.bcd.$bcdist/process_radtags.data.log | cut -d ":" -f 1 )
  head -n $(( $n2 - 2 )) $dmx2.pango${seqs}.bcd.$bcdist/process_radtags.data.log | tail -n +$(( $n1 +1)) >> $data/demultiplex_summary.all.txt
done
head $data/demultiplex_summary.all.txt
wc -l $data/demultiplex_summary.all.txt



#####################################################
#### REDOING DlaA30b	130	TGCTGAA	122	TGTGATAT	9	Pango29 DUE TO TYPO ####
#####################################################

NT=2 ; memo=8
plate=Novaseq_2_2021
savepath=/save/jsalmona/pangolins/RAD/data/
for seqs in $(seq 9) ; do
  seqs=9
  echo "treating rpi DlaA30b pango${seq} second novaseq" 
  fseq=$savepath/211206_A00902_B_L1-2_ASWZ-9-9_R1.fastq.gz
  rseq=$savepath/211206_A00902_B_L1-2_ASWZ-9-9_R2.fastq.gz
  awk -v var="$seq" '$6 == var' $data/RUN2_novaseq_2021.data.info | grep DlaA30b | awk '{print $3 "\t" $1}' > $data/Novaseq_2_2021.pango${seqs}.barcode.Dlaredo 
  bcodefile=$data/Novaseq_2_2021.pango${seqs}.barcode.Dlaredo  ; ddr="no" ; bcdist=1 ; out_dir=$dmx2.pango${seqs}.bcd.${bcdist} # this considers only the first barcode and a barcode distance of 1bp
  #tempbar=$data/tempbarcode #only thing in file is TGCTGAA  DlaA30b
  enzyme1=pstI ; enzyme2=mspI
  mkdir -p $out_dir ; cd $out_dir
  sbatch --cpus-per-task=${NT} --mem=${memo}G  -J pcsrd.$plate.$seq -o $bin/pcsrd.$bcdist.$plate.$seqs.DlaA30b $bin/process_radtags.pg.sh $fseq $rseq $out_dir $bcodefile $ddr $bcdist $enzyme1 $enzyme2
  #$bin/process_radtags.pg.sh $fseq $rseq $out_dir $bcodefile $ddr $bcdist $enzyme1 $enzyme2
  cd $bin
done
squeue -u $user -o "%.15i %.35j %.8u %.2t %.10M %.6D %.6C %.6m %R"



