#!/bin/sh
#!/usr/bin/env Rscript

#called by 01.alignmentstats_coverage.sh

args = commandArgs(trailingOnly=TRUE)

list=args[1]
cov_dir=args[2]
loc=args[3]
setwd(cov_dir)

print(paste("file list is",list))
files<-read.table(paste0(list),header=F)
head(files)


print("creating pdf")
pdf(paste0(loc, ".f1.cov.pdf"),15,15)

par(mfrow=c(3,3))

for(indv in files[,1]){
print(paste("reading individual:",indv))

tab<-read.table(paste0(indv, ".", loc, ".f1.bamhits"),header=F)
tab[,1]<-as.numeric(tab[,1])
head(tab)

tabx<-tab[tab[,1]>1,]
medcovx=median(tabx[,1]) ;

tab2<-tabx[tabx[,1]<medcovx*6,]
maxi=round(quantile(tab2[,1], probs=0.996),digits=1) 
mini=round(quantile(tab2[,1], probs=0.005),digits=1)
medcov=median(tab2[,1])

hist1<-hist(tab2[,1],xlim=c(0,50),ylim=c(0,1000),breaks=seq(1,max(tab2[,1]),1),col="green",xlab="Coverage",main=indv) 

  text(maxi+1,3.5e+3,paste("max cov",maxi),col="blue",srt=90)
  text(15,9e+3,paste(sum(hist1$counts[mini:maxi]),"loci with cov >",mini,"x and <",maxi,"x"),col="darkgreen")
  text(15,1e+4,paste("median cov =",medcov),col="red")
  abline(v=maxi,col="blue")
  text(0,5e+3,paste(length(tab[tab[,1]<2,1]),"loci with cov=1 (not shown)"),col="blue",srt=90)
}
dev.off()

