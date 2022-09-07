#!/bin/sh
#!/usr/bin/env Rscript

#called by coverage.sh which is called by 01.alignmentstats_coverage.sh

args = commandArgs(trailingOnly=TRUE)

indv=args[1]
cov_dir=args[2]

setwd(cov_dir)

tab<-read.table(paste0(indv,".f1.bamhits"),header=F)
tab[,1]<-as.numeric(tab[,1])
head(tab)

if(sum(tab[,1])<4){print("not enough data, quiting"); q(save='no')}

tab2<-tab[tab[,1]>1,]
maxi=round(quantile(tab2[,1], probs=0.990),digits=1) ; print(paste0("quantile 99.0 = ",maxi))
mini=round(quantile(tab2[,1], probs=0.005),digits=1) ; print(paste0("quantile 00.5 = ",mini))
medcov=median(tab2[,1]) ; print(paste0("median = ",medcov))

pdf(paste0(indv,".f1.cov.pdf"))
hist0<-hist(tab2[,1],breaks=seq(1,max(tab2[,1]),1),plot=F)
hist1<-hist(tab2[,1],xlim=c(0,1.2*maxi),ylim=c(0,1.2*max(hist0$counts)),breaks=seq(1,max(tab2[,1]),1),col="green",xlab="Coverage",main=indv) #ylim=c(0,0.05e+3)
  text(maxi+1,0.55*max(hist0$counts),paste("max cov",maxi),col="blue",srt=90)
  text(15,1.05*max(hist0$counts),paste(sum(hist1$counts[mini:maxi]),"loci with cov >",mini,"x and <",maxi,"x"),col="darkgreen")
  text(15,1.1*max(hist0$counts),paste("median cov =",medcov),col="red")#system(paste("echo ",file_map[index]," ",mini," ",maxi," >> ./ind_cov.",gr,".txt",sep=""))
  abline(v=maxi,col="blue")
  #abline(v=mini,col="red")
  text(0,0.55*max(hist0$counts),paste(length(tab[tab[,1]<2,1]),"loci with cov=1 (not shown)"),col="blue",srt=90)
dev.off()

print(paste(sum(hist1$counts[mini:maxi]),"loci with cov >",mini,"x and <",maxi,"x"))
print(paste(length(tab[tab[,1]<2,1]),"loci with cov=1 (not shown)"))
