#!/bin/sh

# module load system/R-3.5.2 libraries/gdal-2.3.0 libraries/proj-4.9.3 system/R-3.5.2 ; R

.libPaths(c( .libPaths(), "/work/pgaubert/softwares/R/3.4.3/lib/", "/work/lchikhi/softwares/R/3.5.2/lib/"))
# install.packages(c('adegenet'), dependencies = TRUE, repos="https://cloud.r-project.org",lib="/work/pgaubert/softwares/R/3.4.3/lib/")

args = commandArgs(trailingOnly=TRUE)
directory<-args[1]
#plink<-args[2]
bindir<-args[2]

library("adegenet")
library("vcfR")

#directory="/work/pgaubert/monk_seals/RAD/results/stacks/X/populations/" ; bindir="/work/pgaubert/monk_seals/RAD/bin"
id="no_filter"

setwd(directory)

# source SNP:error script
source(paste0(bindir,"/SNPs_error.R"))


# Create genlight object based on plink file
#liSNPs<- read.PLINK(file = paste0(directory,"/",plink))

##########################################################
# estimate error rates as a function of depth ############
##########################################################

maxdp=25
err_rates<-as.data.frame(matrix(NA,ncol=(maxdp+1),nrow=5))
for(dp in 1:maxdp){
  vcf <- read.vcfR(paste0(directory,"/populations.snps.dp.",dp,".recode.vcf"), verbose = FALSE)
  liSNPs <- vcfR2genlight(vcf)

  # change names to fit error script requirements
  #liSNPs@ind.names<-gsub("_1","",liSNPs@ind.names)
  liSNPs@ind.names<-gsub("_new","_r",liSNPs@ind.names)

  # Transform param to id (do not include _)
  #id<-gsub(tsv, pattern="(.SNP.SNPs){0,1}(JmonExplo_){0,1}(JmonExplo){0,1}", replacement="")
  id="test1"
  
  # Run funtion to get SNP error rate
  y<-SNP_error(liSNPs=liSNPs, param=id)
  y<-as.data.frame(y)  

  print(dp)  
  print(y)
  err_rates[,(dp+1)]<-y[,2]
}
err_rates[,1]<-y[,1]

write.table(err_rates,file="err_rates",sep="\t",quote=FALSE) 

# plot error rates as a function of depth
pdf("error_rate_x_dp.pdf",6,6)
par(mar=c(5.1,5.1,1.5,1.5))
plot(x=1:maxdp,y=t(err_rates[1,2:(maxdp+1)]),ylim=c(0,0.2),type="l",col=1,lty=1,lwd=1.5,
  xlab="Minimum depth",ylab="Error rate")
for(i in 2:5){
lines(x=1:maxdp,y=t(err_rates[i,2:(maxdp+1)]),col=i,lty=i,lwd=1.5)}
legend("topright",lty=1:5,col=1:5,legend=liSNPs@ind.names[seq(1,9,2)],lwd=1.5)
dev.off()

##########################################################
# estimate error rates as a function of genotype quality #
##########################################################
gqs<-seq(5,40,5)
maxgq<-length(gqs)
err_rates_gq<-as.data.frame(matrix(NA,ncol=(maxgq+1),nrow=5))
for(gq in 1:maxgq){
gqn<-gqs[gq]
vcf <- read.vcfR(paste0(directory,"/populations.snps.gq.",gqn,".recode.vcf"), verbose = FALSE)
liSNPs <- vcfR2genlight(vcf)
#liSNPs@ind.names<-gsub("_1","",liSNPs@ind.names)
liSNPs@ind.names<-gsub("_new","_r",liSNPs@ind.names)
y<-SNP_error(liSNPs=liSNPs, param=id)
y<-as.data.frame(y)  
print(gqn)  
print(y)
err_rates_gq[,(gq+1)]<-y[,2]
}
err_rates_gq[,1]<-y[,1]

write.table(err_rates_gq,file="err_rates_gq",sep="\t",quote=FALSE) 

# plot error rates as a function of genotype quality
pdf("error_rate_x_gq.pdf",6,6)
par(mar=c(5.1,5.1,1.5,1.5))
plot(x=gqs,y=t(err_rates_gq[1,2:(maxgq+1)]),ylim=c(0,0.12),type="l",col=1,lty=1,lwd=1.5,
  xlab="Minimum genotype quality",ylab="Error rate")
for(i in 2:5){
lines(x=gqs,y=t(err_rates_gq[i,2:(maxgq+1)]),col=i,lty=i,lwd=1.5)}
legend("topright",lty=1:5,col=1:5,legend=liSNPs@ind.names[seq(1,9,2)],lwd=1.5)
dev.off()

png("error_rate_all.png",width = 10, height = 5, units = 'in', res = 300)
par(mfrow=c(1,2))

par(mar=c(5.1,5.1,1.5,1.5))
plot(x=1:maxdp,y=t(err_rates[1,2:(maxdp+1)]),ylim=c(0,0.1),type="l",col=1,lty=1,lwd=1,
  xlab="Minimum depth",ylab="Error rate")
for(i in 2:5){
lines(x=1:maxdp,y=t(err_rates[i,2:(maxdp+1)]),col=i,lty=i,lwd=1)}
legend("topright",lty=1:5,col=1:5,legend=liSNPs@ind.names[seq(1,9,2)],lwd=1)

par(mar=c(5.1,5.1,1.5,1.5))
plot(x=gqs,y=t(err_rates_gq[1,2:(maxgq+1)]),ylim=c(0,0.1),type="l",col=1,lty=1,lwd=1,
  xlab="Minimum genotype quality",ylab="Error rate")
for(i in 2:5){
lines(x=gqs,y=t(err_rates_gq[i,2:(maxgq+1)]),col=i,lty=i,lwd=1)}
legend("topright",lty=1:5,col=1:5,legend=liSNPs@ind.names[seq(1,9,2)],lwd=1)
dev.off()






