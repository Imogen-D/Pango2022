

nt=1 ; memo=1 ; srun --cpus-per-task=${nt} --mem=${memo}G --pty bash -i #--x11 



# define paths

user=idumville
#user=jsalmona
res=/work/$user/pangolins/RAD/results
bin=/work/$user/pangolins/RAD/bin
data=/work/$user/pangolins/RAD/data

oldres=/work/jsalmona/pangolins/RAD/results
olddata=/work/jsalmona/pangolins/RAD/data
cd $bin

squeue -u $user -o "%.15i %.30j %.8u %.2t %.10M %.6D %.6C %.6m %R"

resamp_dir=$res/resampling_pa



############################################
####### get metadata information ########## 
############################################
cp $res/stacksautosomes/alllineagesURpopmap.txt $resamp_dir/population_data.txt
pop_data=$resamp_dir/population_data.txt
sed  -i '1i AcessID\tPopulation' $pop_data

############################################
####### Now into R script 06.privatealleleresampling.R ########## 
############################################
module load system/R-3.5.1

