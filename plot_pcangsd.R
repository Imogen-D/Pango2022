#!/bin/sh
#!/usr/bin/env Rscript
# new ### module load system/R-3.5.2 ; R
.libPaths( c( .libPaths(), "/work/lchikhi/softwares/R/3.5.2/lib/","/work/gbesnard/softwares/R/3.5.2/lib/","/work/pgaubert/softwares/R/3.5.2/lib/") )

##### call arguments
args = commandArgs(trailingOnly=TRUE)
loc<-args[2]
prefix_sample_file=args[3]
loci=args[4]
npc=args[5]
# args[1]<-"../../data/genomic_data.info" ; loc="/work/pgaubert/monk_seals/combined_RAD_genomes/results/pcangsd/" ;order_file<-"../../data/genomic_data.info" ; loci="X" ; npc="all" ; prefix_sample_file<-paste0(loci,".mms.combined.all.",npc) ; 

setwd(loc)
#if(loci=="Y"){coln<-7 ; print(loci); print(coln)} 
#if(loci=="X"){coln<-8 ; print(loci);  print(coln)} 
#if(loci=="noX_noMT"){coln<-9 ; print(loci);  print(coln)} 

##for nonscritping

#module load system/R-3.5.2 ; Rscript $bin/plot_pcangsd.R $infofile $pcangsdir $outfile $loc $npc

#infofile<-"/work/idumville/pangolins/RAD/results/ngsadmix.xchrom/sampledataxchrom.txt"
#loc<-"/work/idumville/pangolins/RAD/results/pcangsd.xchrom"
#prefix_sample_file<-"/work/idumville/pangolins/RAD/results/pcangsd.xchrom/rad.ptri.pango.rad.all_xchrom.all"
#loci<-"xchrom"
#npc<-2 #for example


#####call libraries
library(ggplot2); library(grid);library(gridExtra); library(RColorBrewer)

tab<-read.table(paste0(prefix_sample_file,".cov"),header=F)
dim(tab)
tab[1:10,1:10]
m <- as.matrix(tab)
e <- eigen(m)

samples<-read.table(args[1],header=F,sep="\t")#samples<-read.table(infofile,header=F,sep="\t")#; samples<-droplevels(samples[samples[,coln]=="yes",]) # coding for other chromosomes
#samples<-which(samples$V1 == "Y199_1")
ns<-length(samples[,1])
popdat<-samples
summary(popdat)


pdf(paste0(prefix_sample_file,".tests.pdf"),100,100)
plot(1:ns,e$values,pch=8,col="blue",xlab="PCs",ylab="eigenvalues")
abline(h=0,col="red",lty=2)
plot(1:ns,e$values/sum(e$values)*100,pch=8,col="blue",xlab="PCs",ylab="% of variance explained",ylim=c(0,20))

# covariance matricx heatmap
colnames(m)<-paste0(substr(popdat$V2,1,5),"_",popdat$V1) #$V2=pop, $V1=indv
rownames(m)<-paste0(substr(popdat$V2,1,5),"_",popdat$V1)

heatmap(m)
#neighbour joining
#plot(ape::nj(m))
plot(hclust(dist(m), "ave"))
dev.off()


cov_euc_dd<-dist(m)
save(cov_euc_dd,file=paste0(prefix_sample_file,"cov_euc_dist.rda"))

message(as.numeric(popdat$V44))

#pdf(paste0(prefix_sample_file,"pcangsd_first_plots.pdf"),12,12)
#par(mfrow=c(2,2))
#plot(1:ns,e$values/sum(e$values)*100,pch=8,col="blue",xlab="PCs",ylab="% of variance explained",ylim=c(0,70))
#legend("topright",pch=1:length(levels(popdat$V2)),col=1:length(levels(popdat$V2)),legend=levels(popdat$V2),cex=0.5)
#for(i in c(1,3,5,7,9,11,13)){ j=i+1;
#plot(e$vectors[,i:j],xlab=paste0("PC",i,"; ",round(e$values[i]/sum(e$values)*100,2),"%"),ylab=paste0("PC",j,"; ",round(e$values[j]/sum(e$values)*100,2),"%"),lwd=2,col=as.numeric(popdat$V2),pch=as.numeric(popdat$V2))}
#dev.off()

#legend change because only 25 pch vlaues and 31 popns
#pchlist <- sample(1:25, 31, replace=T)


#mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(31)
mycolors <- brewer.pal(6, "Set3")
message(mycolors)
#mypch <- rep(seq.int(0, 25), 2)
mypch <- c(21,22,23,24,25,21)
message(mypch)

pdf(paste0(prefix_sample_file,"pcangsd_first_plots.pdf"),12,12)
par(mfrow=c(2,2))
plot(1:ns,e$values/sum(e$values)*100,pch=8,col="blue",xlab="PCs",ylab="% of variance explained",ylim=c(0,70))
#legend("topright", legend=levels(popdat$V2), cex=0.5, col=c(mycolors), pch=c(mypch))
legend("topright", legend=c("1","2","3", "4", "5", "6"), cex=0.5, col=c(mycolors), pch=c(mypch))

for(i in c(1,3,5,7,9,11,13)){ j=i+1;
plot(e$vectors[,i:j],xlab=paste0("PC",i,"; ", round(e$values[i]/sum(e$values)*100,2),"%"),ylab=paste0("PC",j,"; ",round(e$values[j]/sum(e$values)*100,2),"%"), lwd=0.5, col="black", cex=1.5, bg=c(mycolors)[as.numeric(popdat$V4)], pch=(mypch)[as.numeric(popdat$V4)])} #col=c(mycolors)[as.numeric(popdat$V4)], pch=c(mypch)[as.numeric(popdat$V4)])}
dev.off()
###


#tiff(paste(prefix_sample_file,"pcangsd_first_plots.tiff",sep=""),width = 12, height = 4.5, units = 'in', res = 300,compression ="lzw")
png(paste0(prefix_sample_file,"pcangsd_first_plots.png"),width = 12, height = 4.5, units = 'in', res = 600)#, type="Xlib")
par(mfrow=c(1,3))
par(mar=c(4.1,4.1,3.1,1.1))
#plot(1:ns,e$values/sum(e$values)*100,pch=8,col="blue",xlab="PCs",ylab="% of variance explained",ylim=c(0,15),main="(A)")
#legend("topright", col=c(mycolors), pch=c(mypch),legend=levels(popdat$V2),cex=0.5)
#legend("topright", legend=c("1","2","3", "4", "5", "6"), cex=0.5, col=c(mycolors), pch=c(mypch))
for(i in c(1,3,5)){ j=i+1;
letters<-c("(B)",1,"(C)")
plot(e$vectors[,i:j],xlab=paste0("PC",i,"; ",round(e$values[i]/sum(e$values)*100,2),"%"),ylab=paste0("PC",j,"; ",round(e$values[j]/sum(e$values)*100,2),"%"),lwd=0.5, col="black", cex=1.5, bg=c(mycolors)[as.numeric(popdat$V4)], pch=(mypch)[as.numeric(popdat$V4)])}
dev.off()


#lwd=2,col=c(mycolors)[as.numeric(popdat$V2)], pch=c(mypch)[as.numeric(popdat$V2)],main=letters[i])}
#p1<-ggplot(e$vector, aes(x=V1, y=V2,group=popdat$for.,colour=popdat$for.)) +geom_point() + theme_bw() + theme(legend.position="none")  
#p2<-ggplot(tab, aes(x=V3, y=V4,group=popdat$for.,colour=popdat$for.)) +geom_point() + theme_bw() + theme(legend.position="none")  
#p3<-ggplot(tab, aes(x=V5, y=V6,group=popdat$for.,colour=popdat$for.)) +geom_point() + theme_bw() + theme(legend.position="none")  
#p4<-ggplot(tab, aes(x=V7, y=V8,group=popdat$for.,colour=popdat$for.)) +geom_point() + theme_bw()   
#p4 <- p4 + theme(legend.position = c(0.1, 0.8))
#p <- grid.arrange(p1, p2, p3, p4, nrow = 2)
#print(p)
#
#+ scale_fill_manual(values=1:length(levels(popdat$for.)), 
#                       name="Forest",
#                       breaks=LETTERS[1:(ncol(admix)-2)],
#                       labels=popdat$for.)

quit(save="no")
####Plotting final PCA for all autosomes all lineages

#final pca plot for all lineages

library(ggplot2); library(grid);library(gridExtra); library(RColorBrewer)

loc="/work/idumville/pangolins/RAD/results/pcangsd.autosomes/"
loci="autosomes"
npc="all"
prefix_sample_file="/work/idumville/pangolins/RAD/results/pcangsd.autosomes/rad.ptri.pango.rad.all_autosomes.all"
infofile="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/sampledataautosomes.txt"

setwd(loc)

samples<-read.table(infofile,header=F,sep="\t")
ns<-length(samples[,1])
popdat<-samples
summary(popdat)

popdat$Colour <- "3"
popdat$Colour[which(popdat$V9 %in% c("WAf", "Gha"))] <- "2"
popdat$Colour[which(popdat$V9 == c("WCA"))] <- "6"
popdat$Colour[which(popdat$V9 == c("Gab"))] <- "4"
popdat$Colour[which(popdat$V9 == c("CA"))] <- "4"
popdat$Colour[which(popdat$V2 %in% c("Bioko_Batete", "Douala_Dakat_Market" , "Douala_Central_Market", "Nditam", "Yabassi",
                                                                   "Bayib-Assibong", "Manyemen_Nguti" , "Foumbot", "Bayomen" ))] <- "5"
#west = 2 DG =3 nrthcam=5 south =6 gab and congo =4
urbs <- which(popdat$V10 == "Urban_Market")

popdat <- popdat[-which(samples$V2 == "Banyo"),]
popdat <- popdat[-urbs,]


tab<-read.table(paste0(prefix_sample_file,".cov"),header=F)
tab <- tab[-which(samples$V2 == "Banyo"),-which(samples$V2 == "Banyo")]
tab <- tab[-urbs,-urbs]
dim(tab)
tab[1:10,1:10]
m <- as.matrix(tab)
e <- eigen(m)

mycolors <- brewer.pal(6, "Set3")
message(mycolors)

mypch <- c(21,22,23,24,25,21)
message(mypch)

png(paste0(prefix_sample_file,"Figure2B.png"),width = 12, height = 4, units = 'in', res = 600)#, type="Xlib")
par(mfrow=c(1,3))
par(mar=c(4.1,4.1,3.1,1.1))
#plot(1:ns,e$values/sum(e$values)*100,pch=8,col="blue",xlab="PCs",ylab="% of variance explained",ylim=c(0,15),main="(A)")
#legend("topright", col=c(mycolors), pch=c(mypch),legend=levels(popdat$V2),cex=0.5)
#legend("topright", legend=c("1","2","3", "4", "5", "6"), cex=0.5, col=c(mycolors), pch=c(mypch))
for(i in c(1)){ j=i+1;
plot(e$vectors[,i:j],xlab=paste0("PC",i,"; ",round(e$values[i]/sum(e$values)*100,2),"%"),ylab=paste0("PC",j,"; ",round(e$values[j]/sum(e$values)*100,2),"%"),lwd=0.5, col="black", cex=1.5, bg=c(mycolors)[as.numeric(popdat$Colour)], pch=(mypch)[as.numeric(popdat$Colour)])
mtext("B", 2, adj=1, las=1, line=1, padj=-12, font=2)}
for(i in c(3)){ j=i+1;
plot(e$vectors[,i:j],xlab=paste0("PC",i,"; ",round(e$values[i]/sum(e$values)*100,2),"%"),ylab=paste0("PC",j,"; ",round(e$values[j]/sum(e$values)*100,2),"%"),lwd=0.5, col="black", cex=1.5, bg=c(mycolors)[as.numeric(popdat$Colour)], pch=(mypch)[as.numeric(popdat$Colour)])
mtext("C", 2, adj=1, las=1, line=1, padj=-12, font=2)}      
for(i in c(5)){ j=i+1;
plot(e$vectors[,i:j],xlab=paste0("PC",i,"; ",round(e$values[i]/sum(e$values)*100,2),"%"),ylab=paste0("PC",j,"; ",round(e$values[j]/sum(e$values)*100,2),"%"),lwd=0.5, col="black", cex=1.5, bg=c(mycolors)[as.numeric(popdat$Colour)], pch=(mypch)[as.numeric(popdat$Colour)])
mtext("D", 2, adj=1, las=1, line=1, padj=-12, font=2)
legend("bottomright", legend=c("Western","DG", "Central", "West Cameroon", "South Cameroon"), cex=1, pt.bg=c("#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462"), pch=c(22,23,24,25,21))
}
dev.off()

########################################
##final plot of PCA for WCA samples#####
########################################

library(ggplot2); library(grid);library(gridExtra); library(RColorBrewer)

loc="/work/idumville/pangolins/RAD/results/pcangsd.concat.ruralhighread_wobio"
loci="concat"
npc="ruralhighread"
prefix_sample_file="/work/idumville/pangolins/RAD/results/pcangsd.concat.ruralhighread_wobio/ruralhighread_wobio.pango.rad.autosomes.all"
infofile="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads/sampledataautosomesruralhighread_wobio.txt"

setwd(loc)

samples<-read.table(infofile,header=F,sep="\t")
ns<-length(samples[,1])
popdat<-samples
summary(popdat)
popdat <- popdat[-which(samples$V2 == "Banyo"),]


tab<-read.table(paste0(prefix_sample_file,".cov"),header=F)
tab <- tab[-which(samples$V2 == "Banyo"),-which(samples$V2 == "Banyo")]

dim(tab)
tab[1:10,1:10]
m <- as.matrix(tab)
e <- eigen(m)

mycolors <- brewer.pal(6, "Set3")
message(mycolors)

mypch <- c(21,22,23,24,25,21)
message(mypch)

png(paste0(prefix_sample_file,"Figure3C.png"),width = 4, height = 12, units = 'in', res = 600)#, type="Xlib")
par(mfrow=c(3,1))
par(mar=c(4.1,4.1,3.1,1.1))
#plot(1:ns,e$values/sum(e$values)*100,pch=8,col="blue",xlab="PCs",ylab="% of variance explained",ylim=c(0,15),main="(A)")
#legend("topright", col=c(mycolors), pch=c(mypch),legend=levels(popdat$V2),cex=0.5)
#legend("topright", legend=c("1","2","3", "4", "5", "6"), cex=0.5, col=c(mycolors), pch=c(mypch))
for(i in c(1)){ j=i+1;
plot(e$vectors[,i:j],xlab=paste0("PC",i,"; ",round(e$values[i]/sum(e$values)*100,2),"%"),ylab=paste0("PC",j,"; ",round(e$values[j]/sum(e$values)*100,2),"%"),lwd=0.5, col="black", cex=1.5,   bg=c(mycolors)[as.numeric(popdat$V4)], pch=(mypch)[as.numeric(popdat$V4)])
mtext("C", 2, adj=1, las=1, line=1, padj=-12, font=2)
legend("topleft", legend=c("Gabon","Medouneu","Yaounde", "West Cameroon", "Bioko", "South Cameroon"), cex=1, pt.bg=c(mycolors), pch=c(mypch))}
for(i in c(3)){ j=i+1;
plot(e$vectors[,i:j],xlab=paste0("PC",i,"; ",round(e$values[i]/sum(e$values)*100,2),"%"),ylab=paste0("PC",j,"; ",round(e$values[j]/sum(e$values)*100,2),"%"),lwd=0.5, col="black", cex=1.5, bg=c(mycolors)[as.numeric(popdat$V4)], pch=(mypch)[as.numeric(popdat$V4)])
mtext("D", 2, adj=1, las=1, line=1, padj=-12, font=2)}      
for(i in c(5)){ j=i+1;
plot(e$vectors[,i:j],xlab=paste0("PC",i,"; ",round(e$values[i]/sum(e$values)*100,2),"%"),ylab=paste0("PC",j,"; ",round(e$values[j]/sum(e$values)*100,2),"%"),lwd=0.5, col="black", cex=1.5, bg=c(mycolors)[as.numeric(popdat$V4)], pch=(mypch)[as.numeric(popdat$V4)])
mtext("E", 2, adj=1, las=1, line=1, padj=-12, font=2)} 
dev.off()


