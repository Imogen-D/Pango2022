#!/bin/sh
#!/usr/bin/env Rscript
# srun --mem=3G --pty bash ; module load system/R-3.4.3 ; R

# THIS IS CALLED BY 02b.ploidy.sh

.libPaths( c( .libPaths(), "/work/gbesnard/softwares/R/3.4.3/") )
# library(devtools)
# install.packages("vanquish",repos="https://cloud.r-project.org",lib="/work/gbesnard/softwares/R/3.4.3/")  
# install.packages("readtext",repos="https://cloud.r-project.org",lib="/work/gbesnard/softwares/R/3.4.3/")  
# devtools::install_github("bbanbury/phrynomics", lib="/work/gbesnard/R/3.4.1/lib/")                        
library(vanquish)


args = commandArgs(trailingOnly=TRUE)
vdir<-args[1] ; print(vdir)
rads<-args[2] ; print(rads)

tab<-args[3]
data.info<-read.table(tab,header=F)
X<-args[4]

setwd(vdir)

test.samples<-data.info[,1]  # <-data.info[data.info$ref2=="y",1];

# load 1 sample vcf file
results<-matrix(nrow=length(test.samples),ncol=13)
pdf(paste0("pop.out.Rr0.4.vanquish.summary.depthabv", X, ".pdf"),10,10)
par(mfrow=c(2,2))
 #min coverage considered 
for(i in 1:length(test.samples)){ #5){ #:20){ #
  id=test.samples[i]
  print(id)
  try({
  vcfile=paste(rads,id,".recode.vcf.gz",sep="")
  myvcf<-read_vcf(fn=vcfile, vcffor="GATK")
  myvcf$Content  
  myvcf$VCF[1:20,]
  dim(myvcf$VCF)
  myvcf$file_sample_name<-c(myvcf$Content[10],myvcf$file_sample_name[1])
  myvcf$VCF<-myvcf$VCF[myvcf$VCF$ZG!="Complex" & myvcf$VCF$DP >=X,] # keep only usefull hom and het and discard snp with depth <XX
  result <- defcon(file = myvcf,rmCNV = F)
  print(result$stat); results[i,1:9]<-as.matrix(result$stat[1:9])
  print(result$result); results[i,10:12]<-as.matrix(result$result[1:3])
  tmp <- summary_vcf(vcf = myvcf, ZG = 'het')
  plot(tmp$density, main=paste(id),ylim=c(0,7),xlim=c(0,1))
    abline(v=0.35,col=2,lty=2)
    abline(v=0.65,col=2,lty=2)
    text(0.1,2,paste("X=",X,";",round(100*(sum(tmp$density$y[which(tmp$density$x>=0.35 & tmp$density$x<=0.65)])/sum(tmp$density$y)),1),"%"))
  for(j in 1:4){
    Xs<-c(10,20,30,40)
    myvcf$VCF<-myvcf$VCF[myvcf$VCF$ZG!="Complex" & myvcf$VCF$DP >=Xs[j],] # keep only usefull hom and het and discard snp with depth <XX
    result <- defcon(file = myvcf,rmCNV = F)
    tmp <- summary_vcf(vcf = myvcf, ZG = 'het')
    lines(tmp$density, col=j+1)
    text(0.1,j+2,paste("X=",Xs[j],";",round(100*(sum(tmp$density$y[which(tmp$density$x>=0.35 & tmp$density$x<=0.65)])/sum(tmp$density$y)),1),"%"),col=j+1)
    }
  })
  #results[i,13]<- sum(tmp$density$y[which(tmp$density$x>=0.35 & tmp$density$x<=0.65)])/sum(tmp$density$y)
  #tmpho <- summary_vcf(vcf = myvcf, ZG = 'hom')
  #lines(tmpho$density, col="red")
}
dev.off()


write.table(results,file="vanquish_results.txt",quote=F)

###################################
###########Now just tongues #######
###################################



vanquishdir<-args[1]
rads<-args[2]
tab<-args[3]
data.info<-read.table(tab,header=F)
X<-args[4]

setwd(vanquishdir)

test.samples<-data.info[,1]  # <-data.info[data.info$ref2=="y",1];

# load 1 sample vcf file
results<-matrix(nrow=length(test.samples),ncol=13)
pdf(paste0("pop.out.Rr0.4.vanquish.summary.depthabv.TONGUES", X, ".pdf"),10,10)
par(mfrow=c(2,2))
 #min coverage considered 
for(i in 1:length(test.samples)){ #5){ #:20){ #
  id=test.samples[i]
  print(id)
  try({
  vcfile=paste(rads,id,".recode.vcf.gz",sep="")
  myvcf<-read_vcf(fn=vcfile, vcffor="GATK")
  myvcf$Content  
  myvcf$VCF[1:20,]
  dim(myvcf$VCF)
  myvcf$file_sample_name<-c(myvcf$Content[10],myvcf$file_sample_name[1])
  myvcf$VCF<-myvcf$VCF[myvcf$VCF$ZG!="Complex" & myvcf$VCF$DP >=X,] # keep only usefull hom and het and discard snp with depth <XX
  result <- defcon(file = myvcf,rmCNV = F)
  print(result$stat); results[i,1:9]<-as.matrix(result$stat[1:9])
  print(result$result); results[i,10:12]<-as.matrix(result$result[1:3])
  tmp <- summary_vcf(vcf = myvcf, ZG = 'het')
  plot(tmp$density, main=paste(id),ylim=c(0,7),xlim=c(0,1))
    abline(v=0.35,col=2,lty=2)
    abline(v=0.65,col=2,lty=2)
    text(0.1,2,paste("X=",X,";",round(100*(sum(tmp$density$y[which(tmp$density$x>=0.35 & tmp$density$x<=0.65)])/sum(tmp$density$y)),1),"%"))
  for(j in 1:4){
    Xs<-c(10,20,30,40)
    myvcf$VCF<-myvcf$VCF[myvcf$VCF$ZG!="Complex" & myvcf$VCF$DP >=Xs[j],] # keep only usefull hom and het and discard snp with depth <XX
    result <- defcon(file = myvcf,rmCNV = F)
    tmp <- summary_vcf(vcf = myvcf, ZG = 'het')
    lines(tmp$density, col=j+1)
    text(0.1,j+2,paste("X=",Xs[j],";",round(100*(sum(tmp$density$y[which(tmp$density$x>=0.35 & tmp$density$x<=0.65)])/sum(tmp$density$y)),1),"%"),col=j+1)
    }
  })
  #results[i,13]<- sum(tmp$density$y[which(tmp$density$x>=0.35 & tmp$density$x<=0.65)])/sum(tmp$density$y)
  #tmpho <- summary_vcf(vcf = myvcf, ZG = 'hom')
  #lines(tmpho$density, col="red")
}
dev.off()


write.table(results,file="tongue_vanquish_results.txt",quote=F)


#positive1<-generate_feature(myvcf,mixture=0, hom_p = 0.999, het_p = 0.5, hom_rho = 0.005, het_rho = 0.1, homcut = 0.99, highcut = 0.7, hetcut = 0.3)
#negative1<-generate_feature(vcf_example,mixture=1)
#train_ct(positive1)



