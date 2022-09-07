#!/usr/bin/env Rscript 


args = commandArgs(trailingOnly=TRUE)

mixref <- args[1]
num <- args[2]
suff <- args[3]
cols<-args[4]

#cols=572
#suff="minmac2"
#num=95
#mixref="/work/idumville/pangolins/RAD/results/assignment/rubiasinput95minmac2.txt"


#mixref="rubias100concatfull.txt"
#num=100
#suff="concatfull"

setwd(paste("/work/idumville/pangolins/RAD/results/assignment/", num, suff, "assignments/", sep=""))

message("loading packages")
library("BONE")
library(glmnet)
library("assignPOP")
library("rubias")
library("tidyverse")
library("SNPRelate")
library("radiator")
library(SeqArray)
library(klaR)


Classes = c(rep("character",4),rep("integer",cols)) #this is just to get classes intable from  head -n 1  $as_dir/rubiasmixinputconcatfull.txt | awk --field-separator="\t" "{ print NF }"


message("reading the data")
mixrefData= read.table(mixref,colClasses = Classes,header=T)
MixtureData = mixrefData[which(mixrefData$sample_type == "mixture"),]
ReferenceData =mixrefData[which(mixrefData$sample_type == "reference"),]

#######################################
#BONE
message("BONE")
message("tops of frames")
str(ReferenceData[, 1:6])
str(MixtureData[, 1:6])

message("default imputation")

source("/work/idumville/pangolins/RAD/bin/impute.R") #edited to not remove missing // only if less than two indv (2/169 = 0.0118) impute2
MixtureY = impute2(MixtureData) # Simple marker mode imputation
ReferenceY = impute2(ReferenceData)

MixtureY = MixtureY$Y
ReferenceY = ReferenceY$Y

message("Choose only same set of markers:")
setdiff(rownames(MixtureY),rownames(ReferenceY))

Markers = intersect(rownames(MixtureY),rownames(ReferenceY))

message("trimming frames to common markers")
MixtureY = MixtureY[Markers,]
ReferenceY = ReferenceY[Markers,]

SampleSizeBaselinePop = ncol(ReferenceY)

message("binding frames")
Y = cbind(ReferenceY,MixtureY)

lambda = seq(0.4,0.02,length.out=40)

#Y=unique(Y) #don't want to do this cause more power for snps that are the same

MBapprox = LASSOSolPath(Y,lambda,intercept=T,SampleSizeBaselinePop,Baseline=T)

SampleNames=colnames(Y)

message("solpath approach")
NetworkResultsSolpath = SolPathInference(MBapprox,ReferenceData,SampleNames=SampleNames,
                                         SampleSizeBaselinePop,alpha=0.05)

message("WTA approach")                                       
NetworkResultsWTA = WTAInference(MBapprox,ReferenceData,SampleNames=SampleNames,SampleSizeBaselinePop)


quit(save="no")

message("RUBIAS")
##Rubius
#checking for duplicaes
#  toss them into a function. 

numThreads = RcppParallel::defaultNumThreads()
numThreads

message("checking for matchy pairs")
matchy_pairs <- close_matching_samples(D = mixrefData, #all data 
                                       gen_start_col = 5, 
                                       min_frac_non_miss = 0.85, 
                                       min_frac_matching = 0.94
                                       )

message("extracting duplicate frame")
duplicate <- matchy_pairs %>%
  arrange(desc(num_non_miss), desc(num_match)) #how mnay duplicates 
  
message("saving duplicate frame")
write.table(duplicate, paste(num, suff, "rubias_duplicatesamps.txt", sep = "") , row.name=FALSE,  quote=FALSE)

message("genetic mixture analysis Bayes")
mix_est_Baye <- infer_mixture(reference = ReferenceData, 
                         mixture = MixtureData, 
                         gen_start_col = 5,
                         prelim_reps = 100,
                         prelim_burn_in = 50,
                         method="BR")                # can add priors   


message("saving Bayes mixture analysis frames x3")
write.table(mix_est_Baye$mixing_proportions, paste(num, suff, "rubias_mixingpropBaye.txt", sep = ""), row.name=FALSE,  quote=FALSE)
write.table(mix_est_Baye$indiv_posteriors[,1:9], paste(num, suff, "rubias_indiv_posteriorsBaye.txt", sep = ""), row.name=FALSE,  quote=FALSE)
write.table(mix_est_Baye$mix_prop_traces, paste(num, suff, "rubias_mix_prop_tracesBaye.txt", sep = ""), row.name=FALSE,  quote=FALSE)
write.table(mix_est_Baye$allele_frequencies, paste(num, suff, "rubias_allele_Baye.txt", sep = ""), row.name=FALSE,  quote=FALSE)


message("genetic mixture analysis BS")
mix_est_BS <- infer_mixture(reference = ReferenceData, 
                         mixture = MixtureData, 
                         gen_start_col = 5,
                         method="PB")                # can add priors   


message("saving bootstrap mixture analysis frames x3")
write.table(mix_est_BS$mixing_proportions, paste(num, suff, "rubias_mixingpropBS.txt", sep = ""), row.name=FALSE,  quote=FALSE)
write.table(mix_est_BS$indiv_posteriors[,1:9], paste(num, suff, "rubias_indiv_posteriorsBS.txt", sep = ""), row.name=FALSE,  quote=FALSE)
write.table(mix_est_BS$mix_prop_traces, paste(num, suff, "rubias_mix_prop_tracesBS.txt", sep = ""), row.name=FALSE,  quote=FALSE)
write.table(mix_est_BS$bootstrapped_proportions, paste(num, suff, "rubias_bootstrapprop_BS.txt", sep = ""), row.name=FALSE,  quote=FALSE)

message("assess from whether any of populations")
# get the maximum-a-posteriori population for each individual
map_rows <- mix_est_BS$indiv_posteriors %>%
  group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()  

#message("saving max prior table")
#write.table(map_rows, paste(num, suff, "rubias_maximumprior1.txt", sep = ""))

message("graphing")
#expect normal dist #If you see a z-score to the maximum-a-posteriori population for an individual in your mixture sample that is considerably less than z_scores you saw in the reference, then you might infer that the individual doesn’t actually fit any of the populations in the reference well.
normo <- tibble(z_score = rnorm(1e06))
x <- ggplot(map_rows, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black") 
png(paste(num, suff, "rubias_z_distribution1.png", sep = ""), 750, 500)
x
dev.off()

message("self assingming")
sa_ref <- self_assign(reference = ReferenceData, gen_start_col = 5)

sa_to_repu <- sa_ref %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))

message("saving self assingment table")
write.table(sa_to_repu, paste(num, suff, "rubias_selfassign.txt", sep = ""), row.name=FALSE,  quote=FALSE)

quit(save="no")





##########ASSINGPOP#########

################ AssignPOP ####
training <- paste("/work/idumville/pangolins/RAD/results/assignment/assignPOP/trainingassignpop", num, suff, ".genepop", sep="") #ie training
testing <- paste("/work/idumville/pangolins/RAD/results/assignment/assignPOP/testingassignpop", num, suff, ".genepop", sep="")##ie testing


YourGenepop <- read.Genepop(training, pop.names=c("Sangmelima", "Banyo", "NorthG", "Ekombitie", "GES", "Bayomen", "Yabassi", "SO", "PNCM", "Manengole", "Eseka", "Bioko", "Akonolinga", "Foumbot", "Nditam", "Abongmbang"), haploid = FALSE)
YourGenepopRd <- reduce.allele(YourGenepop, p = 0.95)

#save this so I can come back to it interactively
saveRDS(YourGenepopRd, file = paste(num, suff, "YourGenepopRd.rds", sep=""))

#quit(save="yes")


# Restore the object
YourGenepopRd <- readRDS(file = paste(num, suff, "YourGenepopRd.rds", sep=""))

#clusterEvalQ(cl, .libPaths(c("/home/idumville/R/x86_64-pc-linux-gnu-library/4.1", "/tools/R/R-4.1.2_gcc-9.3.0/lib64/R/library"))


assign.MC( YourGenepopRd, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir=paste(num, suff, "2MC_Result/", sep="") , multiprocess=FALSE)
assign.MC( YourGenepopRd, train.inds=c(1,3,5), train.loci=c(0.1, 0.25, 0.5, 1),
           loci.sample="fst", iterations=30, model="svm", dir=paste(num, suff, "2Num_Indv_MC_Result/", sep=""), multiprocess=FALSE) #, processors = 10) #I don't know what this does whne popn below 3 or 5 indv
assign.kfold( YourGenepopRd, k.fold=c(3, 4, 5), train.loci=c(0.1, 0.25, 0.5, 1), 
             loci.sample="random", model="lda", dir=paste(num, suff, "2kfold_Result/", sep="") , multiprocess=FALSE)  #processors = 10)
             
accuMC1 <- accuracy.MC(dir=paste(num, suff, "2MC_Result/", sep="")) #Use this function for Monte-Carlo cross-validation results

library(ggplot2)
a <- accuracy.plot(accuMC1, pop=c("all","Sangmelima", "Banyo", "NorthG", "Ekombitie", "GES", "Bayomen", "Yabassi", "SO", "PNCM", "Manengole", "Eseka", "Bioko", "Akonolinga", "Foumbot", "Nditam", "Abongmbang")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  #annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci proportions")+ #Add a plot title
  theme(plot.title = element_text(size=20, face="bold")) #Edit plot title text size
ggsave(paste(num, suff, "MC_prop_indv_assignment.pdf", sep=""), a, width = 30,
  height = 10,
  units = c("cm"))
accuMC2 <- accuracy.MC(dir=paste(num, suff, "2Num_Indv_MC_Result/", sep="")) #Use this function for Monte-Carlo cross-validation results
a <- accuracy.plot(accuMC1, pop=c("all","Sangmelima", "Banyo", "NorthG", "Ekombitie", "GES", "Bayomen", "Yabassi", "SO", "PNCM", "Manengole", "Eseka", "Bioko", "Akonolinga", "Foumbot", "Nditam", "Abongmbang")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  #annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle(paste(num, suff, "Monte-Carlo cross-validation using genetic loci number indv", sep=""))+ #Add a plot title
  theme(plot.title = element_text(size=20, face="bold")) #Edit plot title text size
ggsave(paste(num, suff, "MC_num_indv_assignment.pdf", sep=""), a, width = 30,
  height = 10,
  units = c("cm"))
accuKF <- accuracy.kfold(dir = paste(num, suff, "2kfold_Result/", sep="")) #Use this function for K-fold cross-validation results

#defaultparams for now
# 1.Perform assignment test using genetic data and naive Bayes
YourGenepopUnknown <- read.Genepop(testing, haploid = FALSE)
assign.X( x1=YourGenepopRd, x2=YourGenepopUnknown, dir=paste(num, suff,"bayes_AP/", sep=""), model="naiveBayes")
assign.X( x1=YourGenepopRd, x2=YourGenepopUnknown, dir=paste(num, suff,"tree_AP/", sep=""), model="tree")
assign.X( x1=YourGenepopRd, x2=YourGenepopUnknown, dir=paste(num, suff,"randomF_AP/", sep=""), model="randomForest")

# 2.Perform assignment test using integrated data and decision tree
#assign.X( x1=YourIntegrateData, x2=YourIntegrateUnknown, dir=paste(num, suff,"Result-folder3/", sep=""), model="tree")

# 3.Perform assignment test uisng non-genetic data and random forest
#assign.X( x1=morphdf_pop, x2=OtherUnknown, dir="Result-folder5/", model="randomForest")

#give this a go but will probably complain
check.loci(dir = paste(num, suff, "2MC_Result/", sep=""), top.loci = 50)
check.loci(dir = paste(num, suff, "2Num_Indv_MC_Result/", sep=""), top.loci = 50)
membership.plot(dir = paste(num, suff, "2kfold_Result/", sep=""), style = 2)

quit(save="no")

quit(save="no")

message("rubias")
message("genetic mixture analysis")
#mix_est <- infer_mixture(reference = ReferenceData, 
                         mixture = MixtureData, 
                         gen_start_col = 5)                # can add priors   
             
             
message("saving mixture analysis frames x3")
write.table(mix_est$mixing_proportions, paste(num, suff, "rubias_mixingprop1.txt", sep = ""), row.name=FALSE,  quote=FALSE)
write.table(mix_est$indiv_posteriors[,1:9], paste(num, suff, "rubias_indiv_posteriors1.txt", sep = ""), row.name=FALSE,  quote=FALSE)
write.table(mix_est$mix_prop_traces, paste(num, suff, "rubias_mix_prop_traces1.txt", sep = ""), row.name=FALSE,  quote=FALSE)

#for now skipping prior 


# do this from commmand line
check.loci(dir = "MC_Result/", top.loci = 50)
check.loci(dir = "Num_Indv_MC_Result/", top.loci = 50)

membership.plot(dir = "kfold_Result/", style = 2)

#now assignemnt ? aybe from command line as will membership plot for me
YourGenepopUnknown <- read.Genepop( testing )

#defaultparams for now
# 1.Perform assignment test using genetic data and naive Bayes
assign.X( x1=YourGenepopRd, x2=YourGenepopUnknown, dir="Result-folder3/", model="naiveBayes")

# 2.Perform assignment test using integrated data and decision tree
assign.X( x1=YourIntegrateData, x2=YourIntegrateUnknown, dir="Result-folder4/", model="tree")

# 3.Perform assignment test uisng non-genetic data and random forest
assign.X( x1=morphdf_pop, x2=OtherUnknown, dir="Result-folder5/", model="randomForest")



