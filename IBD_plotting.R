#!/bin/sh
#!/usr/bin/env Rscript

#to change to arguments
#from first PCAangsd with all samples



#######################################################
############ Loading Packages  ###################
#######################################################

library(ecodist)
library(sp)
library(vegan)
library(RColorBrewer)
library(ade4)
library(stringr)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

args = commandArgs(trailingOnly=TRUE)

##to change
setwd=args[1]
geodata=args[2]
samplelist=args[3]
DGen_concat=args[4]
DGen_x=args[5]
DGen_auto=args[6]
Cam_samples=args[7]


IBD_dir="/work/idumville/pangolins/RAD/results/IBD/IBD.allsamples" 
setwd(IBD_dir) #output directory
DGen_orig_auto="rad.ptri.pango.rad.all_autosomes.allcov_euc_dist.rda"
DGen_orig_xchrom="rad.ptri.pango.rad.all_xchrom.allcov_euc_dist.rda"

lineage_samples="alllineages.txt"

cam_samples="cameroon.txt"
all_samples="allcameroon.txt"
rural_samples="rural.txt"
geodata="geodata.txt" #coordinate file of all samples


all_names <- readLines(paste0("../samplelists/", all_samples)) #U+R Cameroon
rural_names <- readLines(paste0("../samplelists/", rural_samples)) #Rural Gab and WCA
cam_names <- readLines(paste0("../samplelists/", cam_samples)) #R Cameroon
lineage_names <- readLines(paste0("../samplelists/", lineage_samples)) #Rural all lienages


DGen_orig_auto=load(file = DGen_orig_auto)
DGen_orig_auto=cov_euc_dd
DGen_orig_auto <- as.matrix(DGen_orig_auto)
DGen_orig_xchrom=load(file = DGen_orig_xchrom)
DGen_orig_xchrom=cov_euc_dd
DGen_orig_xchrom <- as.matrix(DGen_orig_xchrom)

#Geographic frames
# geographic distances imported as df
Dgeo_orig <- read.table(paste("../",geodata,sep=""), row.names = 1)

Dgeo_all <- Dgeo_orig[c(which(rownames(Dgeo_orig) %in% all_names)),]
#Dgeo_all_x <- Dgeo_all[-c(which(rownames(Dgeo_all) == "Y199_1")),] #Y199_1 aahhh
Dgeo_rural <- Dgeo_orig[c(which(rownames(Dgeo_orig) %in% rural_names)),]
Dgeo_cam <- Dgeo_orig[c(which(rownames(Dgeo_orig) %in% cam_names)),]
Dgeo_line <- Dgeo_orig[c(which(rownames(Dgeo_orig) %in% lineage_names)),]

Dgeo_all <- as.dist(spDists(as.matrix(Dgeo_all), longlat=T))
#Dgeo_all_x <- as.dist(spDists(as.matrix(Dgeo_all_x), longlat=T))
Dgeo_rural <- as.dist(spDists(as.matrix(Dgeo_rural), longlat=T))
Dgeo_cam <- as.dist(spDists(as.matrix(Dgeo_cam), longlat=T))
Dgeo_line <- as.dist(spDists(as.matrix(Dgeo_line), longlat=T))

#Sorting rownames
DGen_orig_auto_names <- sub("(.*?)(_)", "", x= rownames(DGen_orig_auto))
DGen_orig_auto_names <- str_remove(DGen_orig_auto_names, "Co_") #becasue DR_Co_ kept Co_

DGen_orig_xchrom_names <- sub("(.*?)(_)", "", x= rownames(DGen_orig_xchrom))
DGen_orig_xchrom_names <- str_remove(DGen_orig_xchrom_names, "Co_")#becuase atm concatenated

DGen_all_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% all_names)),c(which(DGen_orig_auto_names %in% all_names))])
DGen_rural_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% rural_names)),c(which(DGen_orig_auto_names %in% rural_names))])
DGen_cam_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% cam_names)),c(which(DGen_orig_auto_names %in% cam_names))])
DGen_line_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% lineage_names)),c(which(DGen_orig_auto_names %in% lineage_names))])

DGen_all_xchrom <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% all_names)),c(which(DGen_orig_xchrom_names %in% all_names))])
DGen_rural_xchrom <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% rural_names)),c(which(DGen_orig_xchrom_names %in% rural_names))])
DGen_cam_xchrom <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% cam_names)),c(which(DGen_orig_xchrom_names %in% cam_names))])
DGen_line_xchrom <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% lineage_names)),c(which(DGen_orig_xchrom_names %in% lineage_names))])




#Estimate the distance at which the IBD is highest ####

source("/work/idumville/pangolins/RAD/bin/R2vsDE_v2.R")
lag=20
minX=50
maxY=0.05

png("IBD.maxdistance.fig1.png",15,10,units="in",res=300)
palette("default")
par(mfrow=c(4,2)) #used to be 1,3 as hist/modern/all
par(mar=c(4.1,4.1,3.1,.1))

mylist<-list(Dgeo=Dgeo_line,Dgen=DGen_line_auto)
Dgen <- DGen_line_auto
Dgeo <- Dgeo_line
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.8,addplot=F,col=1)

mylist<-list(Dgeo=Dgeo_line,Dgen=DGen_line_xchrom)
Dgen <- DGen_line_xchrom
Dgeo <- Dgeo_line
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.8,addplot=F,col=2)

mylist<-list(Dgeo=Dgeo_all,Dgen=DGen_all_auto)
Dgen <- DGen_all_auto
Dgeo <- Dgeo_all
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=maxY,addplot=F,col=3, lty=2)

mylist<-list(Dgeo=Dgeo_all,Dgen=DGen_all_xchrom)
Dgen <- DGen_all_xchrom
Dgeo <- Dgeo_all
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=maxY,addplot=F,col=4, lty=2)

mylist<-list(Dgeo=Dgeo_rural,Dgen=DGen_rural_auto)
Dgen <- DGen_rural_auto
Dgeo <- Dgeo_rural
R2vsDE_2<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=maxY,addplot=F,col=5, lty=3)

mylist<-list(Dgeo=Dgeo_rural,Dgen=DGen_rural_xchrom)
Dgen <- DGen_rural_xchrom
Dgeo <- Dgeo_rural
R2vsDE_2<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=maxY,addplot=F,col=6, lty=3)

mylist<-list(Dgeo=Dgeo_cam,Dgen=DGen_cam_auto)
Dgen <- DGen_cam_auto
Dgeo <- Dgeo_cam
R2vsDE_3<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=maxY,addplot=F,col=1, lty=4)

mylist<-list(Dgeo=Dgeo_cam,Dgen=DGen_cam_xchrom)
Dgen <- DGen_cam_xchrom
Dgeo <- Dgeo_cam
R2vsDE_3<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=maxY,addplot=F,col=2, lty=4)

legend("topleft",lty=c(1,1,2,2,3,3,4,4),col=c(1:6,1,2),legend=c("All Autosomes", "All X", "Cameroon R & U Autosomes", "Cameroon R & U X", "WCA & Gab Rural Autosomes", "WCA & Gab Rural X", "Cameroon Rural Autosomes", "Cameroon Rural X"),title="Gen.Data",cex=0.8)

dev.off()

#Highest distance = largest R2


#Mantel correlograms

#making own classes ; change dgeo and dgen for different populations
#cameroon
mc.all.auto<-mantel.correlog(D.eco=DGen_all_auto,D.geo=Dgeo_all,n.class = 50,nperm = 999)
mc.rural.auto<-mantel.correlog(D.eco=DGen_rural_auto,D.geo =Dgeo_rural,n.class = 50,nperm = 999)
mc.cam.auto<-mantel.correlog(D.eco=DGen_cam_auto,D.geo =Dgeo_cam,n.class = 50,nperm = 999)
mc.all.x<-mantel.correlog(D.eco=DGen_all_xchrom,D.geo =Dgeo_all,n.class = 50,nperm = 999)
mc.rural.x<-mantel.correlog(D.eco=DGen_rural_xchrom,D.geo =Dgeo_rural,n.class = 50,nperm = 999)
mc.cam.x<-mantel.correlog(D.eco=DGen_cam_xchrom,D.geo =Dgeo_cam,n.class = 50,nperm = 999)
mc.line.auto<-mantel.correlog(D.eco=DGen_line_auto,D.geo=Dgeo_line,n.class = 50,nperm = 999)
mc.line.x<-mantel.correlog(D.eco=DGen_line_xchrom,D.geo=Dgeo_line,n.class = 50,nperm = 999)

png("IBD.mantel_corr.fig2.png",20,16,units="in",res=300)
palette("default")
par(mfrow=c(4,2))
par(mar=c(4.1,4.1,.1,.1))
#plot(mc.all,ylim=c(-1,1))
#plot All Lineages
plot(mc.line.auto$mantel.res[,1],mc.line.auto$mantel.res[,3],ylim=c(-.25,+.4),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,2500));
axis(1, las=2)
mtext("  A) Autosome Rural All Lineages",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.line.auto$mantel.res[,1],mc.line.auto$mantel.res[,3])
for(i in 1:length(mc.line.auto$mantel.res[,1])){
  if(is.na(mc.line.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.line.auto$mantel.res[i,5]<.05){points(mc.line.auto$mantel.res[i,1],mc.line.auto$mantel.res[i,3],pch=22,col=1,bg=1)}}}
plot(mc.line.x$mantel.res[,1],mc.line.x$mantel.res[,3],ylim=c(-.25,+.4),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,2500));
axis(1, las=2)
mtext("  B) X Rural All Lineages",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.line.x$mantel.res[,1],mc.line.x$mantel.res[,3])
for(i in 1:length(mc.line.x$mantel.res[,1])){
  if(is.na(mc.line.x$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.line.x$mantel.res[i,5]<.05){points(mc.line.x$mantel.res[i,1],mc.line.x$mantel.res[i,3],pch=22,col=1,bg=1)}}}
#plot Rural+Urban
plot(mc.all.auto$mantel.res[,1],mc.all.auto$mantel.res[,3],ylim=c(-.1,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,500));
axis(1, las=2)
mtext("  C) Autosome rural & urban Cameroon",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.all.auto$mantel.res[,1],mc.all.auto$mantel.res[,3])
for(i in 1:length(mc.all.auto$mantel.res[,1])){
  if(is.na(mc.all.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.all.auto$mantel.res[i,5]<.05){points(mc.all.auto$mantel.res[i,1],mc.all.auto$mantel.res[i,3],pch=22,col=1,bg=1)}}}
plot(mc.all.x$mantel.res[,1],mc.all.x$mantel.res[,3],ylim=c(-.1,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,500));
axis(1, las=2)
mtext("  D) X rural & urban Cameroon",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.all.x$mantel.res[,1],mc.all.x$mantel.res[,3])
for(i in 1:length(mc.all.x$mantel.res[,1])){
  if(is.na(mc.all.x$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.all.x$mantel.res[i,5]<.05){points(mc.all.x$mantel.res[i,1],mc.all.x$mantel.res[i,3],pch=22,col=1,bg=1)}}}
#plot Rural
plot(mc.rural.auto$mantel.res[,1],mc.rural.auto$mantel.res[,3],ylim=c(-.1,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,500));
axis(1, las=2)
mtext("  E) Autosome rural",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.rural.auto$mantel.res[,1],mc.rural.auto$mantel.res[,3])
for(i in 1:length(mc.rural.auto$mantel.res[,1])){
  if(is.na(mc.rural.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.rural.auto$mantel.res[i,5]<.05){points(mc.rural.auto$mantel.res[i,1],mc.rural.auto$mantel.res[i,3],pch=22,col=1,bg=1)}}}
plot(mc.rural.x$mantel.res[,1],mc.rural.x$mantel.res[,3],ylim=c(-.1,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,500));
axis(1, las=2)
mtext("  E) X rural",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.rural.x$mantel.res[,1],mc.rural.x$mantel.res[,3])
for(i in 1:length(mc.rural.x$mantel.res[,1])){
  if(is.na(mc.rural.x$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.rural.x$mantel.res[i,5]<.05){points(mc.rural.x$mantel.res[i,1],mc.rural.x$mantel.res[i,3],pch=22,col=1,bg=1)}}}
# plot cameroon
plot(mc.cam.auto$mantel.res[,1],mc.cam.auto$mantel.res[,3],ylim=c(-.1,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,500));
axis(1, las=2)
mtext("  F) Autosome Cameroon",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.cam.auto$mantel.res[,1],mc.cam.auto$mantel.res[,3])
for(i in 1:length(mc.cam.auto$mantel.res[,1])){
  if(is.na(mc.cam.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.cam.auto$mantel.res[i,5]<.05){points(mc.cam.auto$mantel.res[i,1],mc.cam.auto$mantel.res[i,3],pch=22,col=1,bg=1)}}}
plot(mc.cam.x$mantel.res[,1],mc.cam.x$mantel.res[,3],ylim=c(-.1,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,500));
axis(1, las=2)
mtext("  G) X Cameroon ",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.cam.x$mantel.res[,1],mc.cam.x$mantel.res[,3])
for(i in 1:length(mc.cam.x$mantel.res[,1])){
  if(is.na(mc.cam.x$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.cam.x$mantel.res[i,5]<.05){points(mc.cam.x$mantel.res[i,1],mc.cam.x$mantel.res[i,3],pch=22,col=1,bg=1)}}}
  
dev.off()

#Plot mantel correlograms with concat, auto, x and cameroon overlaid

png("IBD.mantel_corr.XACam.fig3.png",10,4,units="in",res=300)
palette(brewer.pal(n = 8, name = "Set2"))
par(mfrow=c(1,1))
par(mar=c(4.1,4.1,.1,.1))
#plot(mc.all,ylim=c(-1,1))

#plot rural autosomes
plot(mc.rural.auto$mantel.res[,1],mc.rural.auto$mantel.res[,3],ylim=c(-.25,+.65),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=3,xaxt="n",xlim=c(0,400));
axis(1, las=2)
abline(h=0,col="red")
lines(mc.rural.auto$mantel.res[,1],mc.rural.auto$mantel.res[,3], col=1, lty=2)
for(i in 1:length(mc.rural.auto$mantel.res[,1])){
  if(is.na(mc.rural.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.rural.auto$mantel.res[i,5]<.05){points(mc.rural.auto$mantel.res[i,1],mc.rural.auto$mantel.res[i,3],pch=22,col=1,bg=1)}}}

#plot rural x
points(mc.rural.x$mantel.res[,1],mc.rural.x$mantel.res[,3],ylim=c(-.25,+.65),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=23,col=3,xaxt="n",xlim=c(0,3000));
axis(1, las=2)
lines(mc.rural.x$mantel.res[,1],mc.rural.x$mantel.res[,3], col=2, lty=3)
for(i in 1:length(mc.rural.x$mantel.res[,1])){
  if(is.na(mc.rural.x$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.rural.x$mantel.res[i,5]<.05){points(mc.rural.x$mantel.res[i,1],mc.rural.x$mantel.res[i,3],pch=22,col=2,bg=2)}}}  
 
legend("bottomright",legend=c("WCA & Gab Rural Autosomes", "WCA & Gab Rural X"),lty=c(2, 3),pch=c(22, 23),col=c(1,2),bty="n")
dev.off()

#######Estimate IBD

ibd.auto.all <- mantel.randtest(DGen_all_auto,Dgeo_all,nrepet = 9999)
ibd.auto.rural <- mantel.randtest(DGen_rural_auto,Dgeo_rural,nrepet = 9999)
ibd.auto.cam <- mantel.randtest(DGen_cam_auto,Dgeo_cam,nrepet = 9999)
ibd.auto.lineages <- mantel.randtest(DGen_line_auto,Dgeo_line,nrepet = 9999)

#### Final figure IBD
#plotting autosome urban, rural and cameroon
png("IBD.main_fig_mantel_corr.fig4_CameronRU.png",16,8,units="in",res=300)
palette("default")
palette(brewer.pal(n = 8, name = "Set2"))
cols <- add.alpha(brewer.pal(n = 8, name = "Set2"), 0.5)

par(mfrow=c(1,1))
par(mar=c(4.1,5.1,.1,.1))

# IBD corr + hist
X<-0.9*max(Dgeo_all)
X2<-0.65*max(Dgeo_all)
plot(Dgeo_all, DGen_all_auto, pch=18,cex=.6, col=cols[1], xlab ="Geographic distance (Km)", #cols[3]
    ylab="Genetic distance",main="",xaxt="n")
axis(side = 1,las=3,cex=0.9)
abline(lm(DGen_all_auto~Dgeo_all), lty=2,col=cols[1]) #cols[3]

points(Dgeo_cam, DGen_cam_auto, pch=18, cex=.6, col=2)
abline(lm(DGen_cam_auto~Dgeo_cam), lty=3,col=cols[2])

text(X,0.2*max(DGen_cam_auto),"R & U", cex=.8)
text(X,0.1*max(DGen_cam_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.auto.all$obs)^2,3))),cex=.8,col="black")
text(X,0.0*max(DGen_cam_auto),substitute(italic(p)<=rr,list(rr=round(ibd.auto.all$pvalue,4))),cex=.8,col="black")

text(X2,0.2*max(DGen_cam_auto),"Rural", cex=.8)
text(X2,0.1*max(DGen_cam_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.auto.cam$obs)^2,3))),cex=.8,col="black")
text(X2,0.0*max(DGen_cam_auto),substitute(italic(p)<=rr,list(rr=round(ibd.auto.cam$pvalue,4))),cex=.8,col="black")

legend("topright",pch=c(18,18),col=c(1,2),legend=c("R & U Cameroon Autosome", "R Cameroon Autosome"),title="Dataset",cex=0.8)

dev.off()


png("IBD.main_fig_mantel_corr.fig4_AllLineages.png",16,8,units="in",res=300)
palette("default")
palette(brewer.pal(n = 8, name = "Set2"))
cols <- add.alpha(brewer.pal(n = 8, name = "Set2"), 0.5)

par(mfrow=c(1,1))
par(mar=c(4.1,5.1,.1,.1))

# IBD with x3 datasets
X1<-0.85*max(Dgeo_line)
X<-0.75*max(Dgeo_line)
X2<-0.65*max(Dgeo_line)
plot(Dgeo_line, DGen_line_auto, pch=18,cex=.6, col=cols[3], xlab ="Geographic distance (Km)", #cols[3]
    ylab="Genetic distance",main="",xaxt="n")
axis(side = 1,las=3,cex=0.9)
abline(lm(DGen_line_auto~Dgeo_line), lty=2,col=cols[3]) #cols[3]

points(Dgeo_rural, DGen_rural_auto, pch=18, cex=.6, col=1)
abline(lm(DGen_rural_auto~Dgeo_rural), lty=3,col=cols[1])

points(Dgeo_cam, DGen_cam_auto, pch=18, cex=.6, col=2)
abline(lm(DGen_cam_auto~Dgeo_cam), lty=3,col=cols[2])

text(X1,0.2*max(DGen_line_auto),"All Lineages Rural", cex=.8)
text(X1,0.1*max(DGen_line_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.auto.lineages$obs)^2,3))),cex=.8,col="black")
text(X1,0.0*max(DGen_line_auto),substitute(italic(p)<=rr,list(rr=round(ibd.auto.lineages$pvalue,4))),cex=.8,col="black")

text(X,0.2*max(DGen_line_auto),"Gabon & Cameroon Rural", cex=.8)
text(X,0.1*max(DGen_line_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.auto.rural$obs)^2,3))),cex=.8,col="black")
text(X,0.0*max(DGen_line_auto),substitute(italic(p)<=rr,list(rr=round(ibd.auto.rural$pvalue,4))),cex=.8,col="black")

text(X2,0.2*max(DGen_line_auto),"Cameroon Rural", cex=.8)
text(X2,0.1*max(DGen_line_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.auto.cam$obs)^2,3))),cex=.8,col="black")
text(X2,0.0*max(DGen_line_auto),substitute(italic(p)<=rr,list(rr=round(ibd.auto.cam$pvalue,4))),cex=.8,col="black")

legend("topright",pch=c(18,18),col=c(3,1,2),legend=c("All Lineages Rural Autosome", "WCA & Gab Rural Autosome", "Cameroon Rural Autosome"),title="Dataset",cex=0.8)

dev.off()


#######################################################
############ Extracting indivudals above thresholds ###################
#######################################################


DGen_temp <- as.matrix(DGen_rural_auto)
DGen_temp[upper.tri(DGen_temp)] <- NA
library(reshape2)
df <- melt(as.matrix(DGen_temp), varnames = c("row", "col"))
df <- df[which(df$value > 4),]
distdf <- dcast(df, row~col)
write.table(distdf, file = "./Gab_Cam_BiggerDB.txt", sep="\t")


DGen_temp <- as.matrix(DGen_all_auto)
DGen_temp[upper.tri(DGen_temp)] <- NA
DGen_temp[DGen_temp < 1.5] <- NA
DGen_temp <-DGen_temp[,colSums(is.na(DGen_temp)) != nrow(DGen_temp)]
DGen_temp <-DGen_temp[rowSums(is.na(DGen_temp)) != ncol(DGen_temp), ]
write.table(DGen_temp, file = "./All_extracted_biggerGB.txt", sep="\t")
DGen_all_auto

DGen_concat_cameroon <-DGen_concat_cameroon[colSums(is.na(DGen_concat_cameroon )) != nrow(DGen_concat_cameroon ), ]

write.table(DGen_concat_cameroon, file = "./Cameroon_extracted_biggerGD.txt", sep="\t")


#######################################################
############ Lineage Plotting ###################
#######################################################

IBD_dir="/work/idumville/pangolins/RAD/results/IBD/IBD.allsamples" 
setwd(IBD_dir) #output directory
DGen_orig_auto="rad.ptri.pango.rad.all_autosomes.allcov_euc_dist.rda"
DGen_orig_xchrom="rad.ptri.pango.rad.all_xchrom.allcov_euc_dist.rda"


WCA_samples="WCA.txt"
bio_samples="Bioko.txt"
CA_samples="CA.txt"
DG_samples="DG.txt"
Gab_samples="Gab.txt"
WCAGab_samples="WCAGab.txt"
WAf_samples="WAf.txt"
WCABio_samples="WCABioko.txt"
lineage_samples="alllineages.txt"


geodata="geodata.txt" #coordinate file of all samples


DGen_orig_auto=load(file = DGen_orig_auto)
DGen_orig_auto=cov_euc_dd
DGen_orig_auto <- as.matrix(DGen_orig_auto)
DGen_orig_xchrom=load(file = DGen_orig_xchrom)
DGen_orig_xchrom=cov_euc_dd
DGen_orig_xchrom <- as.matrix(DGen_orig_xchrom)

#Sorting rownames because concatenated
DGen_orig_auto_names <- sub("(.*?)(_)", "", x= rownames(DGen_orig_auto))
DGen_orig_auto_names <- str_remove(DGen_orig_auto_names, "Co_") #becasue DR_Co_ kept Co_
DGen_orig_xchrom_names <- sub("(.*?)(_)", "", x= rownames(DGen_orig_xchrom))
DGen_orig_xchrom_names <- str_remove(DGen_orig_xchrom_names, "Co_")

#Geographic frames
# geographic distances imported as df
Dgeo_orig <- read.table(paste("../samplelists/",geodata,sep=""), row.names = 1)

Dgeo_WCA <- as.dist(spDists(as.matrix(Dgeo_orig[c(which(rownames(Dgeo_orig) %in% readLines(paste0("../samplelists/", WCA_samples)))),]), longlat=T))
Dgeo_bio <- as.dist(spDists(as.matrix(Dgeo_orig[c(which(rownames(Dgeo_orig) %in% readLines(paste0("../samplelists/", bio_samples)))),]), longlat=T))
Dgeo_CA <- as.dist(spDists(as.matrix(Dgeo_orig[c(which(rownames(Dgeo_orig) %in% readLines(paste0("../samplelists/", CA_samples)))),]), longlat=T))
Dgeo_DG <- as.dist(spDists(as.matrix(Dgeo_orig[c(which(rownames(Dgeo_orig) %in% readLines(paste0("../samplelists/", DG_samples)))),]), longlat=T))
Dgeo_Gab <- as.dist(spDists(as.matrix(Dgeo_orig[c(which(rownames(Dgeo_orig) %in% readLines(paste0("../samplelists/", Gab_samples)))),]), longlat=T))
Dgeo_WCAGab <- as.dist(spDists(as.matrix(Dgeo_orig[c(which(rownames(Dgeo_orig) %in% readLines(paste0("../samplelists/", WCAGab_samples)))),]), longlat=T))
Dgeo_WCABio <- as.dist(spDists(as.matrix(Dgeo_orig[c(which(rownames(Dgeo_orig) %in% readLines(paste0("../samplelists/", WCABio_samples)))),]), longlat=T))
Dgeo_WAf <- as.dist(spDists(as.matrix(Dgeo_orig[c(which(rownames(Dgeo_orig) %in% readLines(paste0("../samplelists/", WAf_samples)))),]), longlat=T))
Dgeo_all <- as.dist(spDists(as.matrix(Dgeo_orig[c(which(rownames(Dgeo_orig) %in% readLines(paste0("../samplelists/", lineage_samples)))),]), longlat=T))

#For now just pure lineages
DGen_WCA_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", WCA_samples)))),c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", WCA_samples))))])
DGen_WAf_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", WAf_samples)))),c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", WAf_samples))))])
DGen_Gab_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", Gab_samples)))),c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", Gab_samples))))])
DGen_CA_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", CA_samples)))),c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", CA_samples))))])
DGen_DG_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", DG_samples)))),c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", DG_samples))))])
DGen_all_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", lineage_samples)))),c(which(DGen_orig_auto_names %in% readLines(paste0("../samplelists/", lineage_samples))))])


DGen_WCA_x <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", WCA_samples)))),c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", WCA_samples))))])
DGen_WAf_x <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", WAf_samples)))),c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", WAf_samples))))])
DGen_Gab_x <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", Gab_samples)))),c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", Gab_samples))))])
DGen_CA_x <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", CA_samples)))),c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", CA_samples))))])
DGen_DG_x <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", DG_samples)))),c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", DG_samples))))])
DGen_all_x <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", lineage_samples)))),c(which(DGen_orig_xchrom_names %in% readLines(paste0("../samplelists/", lineage_samples))))])

#Estimate the distance at which the IBD is highest ####

source("/work/idumville/pangolins/RAD/bin/R2vsDE_v2.R")
lag=20
minX=50
maxY=0.05

#WCA/CA/DG/Gab/WAf/all

png("IBD.maxdistance.seplineages.fig1.png",15,10,units="in",res=300)
palette("default")
palette(brewer.pal(n = 8, name = "Set2"))
par(mfrow=c(1,2))
par(mar=c(4.1,4.1,3.1,.1))

#mylist<-list(Dgeo=Dgeo_all,Dgen=DGen_all_auto)
#Dgen <- DGen_all_auto
#Dgeo <- Dgeo_all
#R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.8,addplot=F,col=1)

mylist<-list(Dgeo=Dgeo_CA,Dgen=DGen_CA_auto)
Dgen <- DGen_CA_auto
Dgeo <- Dgeo_CA
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.6,addplot=F,col=2)

mylist<-list(Dgeo=Dgeo_WAf,Dgen=DGen_WAf_auto)
Dgen <- DGen_WAf_auto
Dgeo <- Dgeo_WAf
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.6,addplot=T,col=3)

mylist<-list(Dgeo=Dgeo_WCA,Dgen=DGen_WCA_auto)
Dgen <- DGen_WCA_auto
Dgeo <- Dgeo_WCA
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.6,addplot=T,col=4)


mylist<-list(Dgeo=Dgeo_DG,Dgen=DGen_DG_auto)
Dgen <- DGen_DG_auto
Dgeo <- Dgeo_DG
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.6,addplot=T,col=5)

mylist<-list(Dgeo=Dgeo_Gab,Dgen=DGen_Gab_auto)
Dgen <- DGen_Gab_auto
Dgeo <- Dgeo_Gab
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.6,addplot=T,col=6)

legend("topright",lty=1,col=c(2:6),legend=c("CA", "WAf", "WCA", "DG",  "Gab"),title="Autosome Dataset",cex=0.8)

#Now X
#mylist<-list(Dgeo=Dgeo_all,Dgen=DGen_all_x)
#Dgen <- DGen_all_x
#Dgeo <- Dgeo_all
#R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.8,addplot=F,col=1)

mylist<-list(Dgeo=Dgeo_CA,Dgen=DGen_CA_x)
Dgen <- DGen_CA_x
Dgeo <- Dgeo_CA
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.6,addplot=F,col=2)

mylist<-list(Dgeo=Dgeo_WAf,Dgen=DGen_WAf_x)
Dgen <- DGen_WAf_x
Dgeo <- Dgeo_WAf
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.6,addplot=T,col=3)

mylist<-list(Dgeo=Dgeo_WCA,Dgen=DGen_WCA_x)
Dgen <- DGen_WCA_x
Dgeo <- Dgeo_WCA
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.6,addplot=T,col=4)


mylist<-list(Dgeo=Dgeo_DG,Dgen=DGen_DG_x)
Dgen <- DGen_DG_x
Dgeo <- Dgeo_DG
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.6,addplot=T,col=5)

mylist<-list(Dgeo=Dgeo_Gab,Dgen=DGen_Gab_x)
Dgen <- DGen_Gab_x
Dgeo <- Dgeo_Gab
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.6,addplot=T,col=6)

legend("topright",lty=1,col=c(2:6),legend=c("CA", "WAf","WCA", "DG", "Gab"),title="X Chromosome Dataset",cex=0.8)

dev.off()

#######Estimate IBD

ibd.auto.all <- mantel.randtest(DGen_all_auto,Dgeo_all,nrepet = 9999)
ibd.auto.WCA <- mantel.randtest(DGen_WCA_auto,Dgeo_WCA,nrepet = 9999)
ibd.auto.CA <- mantel.randtest(DGen_CA_auto,Dgeo_CA,nrepet = 9999)
ibd.auto.WAf <- mantel.randtest(DGen_WAf_auto,Dgeo_WAf,nrepet = 9999)
ibd.auto.Gab <- mantel.randtest(DGen_Gab_auto,Dgeo_Gab,nrepet = 9999)
ibd.auto.DG <- mantel.randtest(DGen_DG_auto,Dgeo_DG,nrepet = 9999)

#### Final figure IBD
#plotting autosome urban, rural and cameroon
png("IBD.main_fig_mantel_corr.seplineages.fig4.png",10,5,units="in",res=300)
palette("default")
palette(brewer.pal(n = 8, name = "Set3"))
cols <- add.alpha(brewer.pal(n = 8, name = "Set3"), 0.2)

par(mfrow=c(1,1))
par(mar=c(4.1,5.1,.1,.1))

# IBD corr + hist
X<-0.9*max(Dgeo_WCA)
X2<-0.7*max(Dgeo_WCA)
X3<-0.5*max(Dgeo_WCA)
X4<-0.3*max(Dgeo_WCA)
X5<-0.1*max(Dgeo_WCA)
plot(Dgeo_WCA, DGen_WCA_auto,pch=21,cex=.6, col="grey20", lwd=0.2, bg=cols[5], xlab ="Geographic distance (Km)",ylab="Genetic distance",main="",xaxt="n", ylim=c(0,3),xlim=c(0,1000))
#plot(Dgeo_all, DGen_all_auto, pch=18,cex=.6, col=cols[1], xlab ="Geographic distance (Km)", #cols[3]
#    ylab="Genetic distance",main="",xaxt="n")
axis(side = 1,las=3,cex=0.9)
#abline(lm(DGen_all_auto~Dgeo_all), lty=2,col=cols[1]) #cols[3]

#points(Dgeo_WCA, DGen_WCA_auto, pch=18, cex=.6, col=cols[4])
abline(lm(DGen_WCA_auto~Dgeo_WCA), lty=3,col=5)
points(Dgeo_CA, DGen_CA_auto, pch=21, cex=.6, bg=4, col="grey20", lwd=0.2)
abline(lm(DGen_CA_auto~Dgeo_CA), lty=3,col=4)
points(Dgeo_WAf, DGen_WAf_auto, pch=21, cex=.6, col="grey20", bg=2, lwd=0.2)
abline(lm(DGen_WAf_auto~Dgeo_WAf), lty=3,col=2)
points(Dgeo_DG, DGen_DG_auto, pch=21, cex=.6, col="grey20",bg=3, lwd=0.2)
abline(lm(DGen_DG_auto~Dgeo_DG), lty=3,col=3)
points(Dgeo_Gab, DGen_Gab_auto, pch=21, cex=.6, col="grey20", bg=7, lwd=0.2)
abline(lm(DGen_Gab_auto~Dgeo_Gab), lty=3,col=7)


text(X,0.2*max(DGen_WCA_auto),"CA", cex=.6)
text(X,0.1*max(DGen_WCA_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.auto.CA$obs)^2,3))),cex=.6,col="black")
text(X,0.0*max(DGen_WCA_auto),substitute(italic(p)<=rr,list(rr=round(ibd.auto.CA$pvalue,4))),cex=.6,col="black")

text(X2,0.2*max(DGen_WCA_auto),"WAf", cex=.6)
text(X2,0.1*max(DGen_WCA_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.auto.WAf$obs)^2,3))),cex=.6,col="black")
text(X2,0.0*max(DGen_WCA_auto),substitute(italic(p)<=rr,list(rr=round(ibd.auto.WAf$pvalue,4))),cex=.6,col="black")

text(X3,0.2*max(DGen_WCA_auto),"WCA", cex=.6)
text(X3,0.1*max(DGen_WCA_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.auto.WCA$obs)^2,3))),cex=.6,col="black")
text(X3,0.0*max(DGen_WCA_auto),substitute(italic(p)<=rr,list(rr=round(ibd.auto.WCA$pvalue,4))),cex=.6,col="black")

text(X4,0.2*max(DGen_WCA_auto),"DG", cex=.6)
text(X4,0.1*max(DGen_WCA_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.auto.DG$obs)^2,3))),cex=.6,col="black")
text(X4,0.0*max(DGen_WCA_auto),substitute(italic(p)<=rr,list(rr=round(ibd.auto.DG$pvalue,4))),cex=.6,col="black")

legend("topright",pch=21, pt.bg=c(4,2,5,3,7), col="grey20", legend=c("CA", "WAf","WCA", "DG", "Gab"),title="Lineage",cex=0.8)

dev.off()

#### Supp Figure 12 only WCA IBD

png("Supp_FigS12_IBD.main_fig_mantel_corr.WCA.png",10,5,units="in",res=300)
palette("default")
palette(brewer.pal(n = 8, name = "Set3"))
cols <- add.alpha(brewer.pal(n = 8, name = "Set3"), 0.2)

par(mfrow=c(1,1))
par(mar=c(4.1,5.1,.1,.1))

X<-0.9*max(Dgeo_WCA)

plot(Dgeo_WCA, DGen_WCA_auto,pch=21,cex=.6, col="grey20", lwd=0.2, bg=cols[5], xlab ="Geographic distance (Km)",ylab="Genetic distance",main="",xaxt="n", ylim=c(0,1.5),xlim=c(0,750))
axis(side = 1,las=3,cex=0.9)
abline(lm(DGen_WCA_auto~Dgeo_WCA), lty=3,col=5)

text(X,0.1*max(DGen_WCA_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.auto.WCA$obs)^2,3))),cex=.6,col="black")
text(X,0.0*max(DGen_WCA_auto),substitute(italic(p)<=rr,list(rr=round(ibd.auto.WCA$pvalue,4))),cex=.6,col="black")

#legend("topright",pch=21, pt.bg=c(4,2,5,3,7), col="grey20", legend=c("CA", "WAf","WCA", "DG", "Gab"),title="Lineage",cex=0.8)

dev.off()


### Ecodist Plot #######
#No Gab because only 2 smaples
stepsize=20
mcc.all<-mgram(species.d=DGen_all_auto, space.d=Dgeo_all, stepsize=stepsize)# nclass=25)
mcc.WCA<-mgram(species.d=DGen_WCA_auto, space.d=Dgeo_WCA, stepsize=stepsize)# nclass=25)
mcc.CA<-mgram(species.d=DGen_CA_auto, space.d=Dgeo_CA, stepsize=stepsize)# nclass=25)
mcc.WAf<-mgram(species.d=DGen_WAf_auto, space.d=Dgeo_WAf, stepsize=stepsize)# nclass=25)
mcc.DG<-mgram(species.d=DGen_DG_auto, space.d=Dgeo_DG, stepsize=stepsize)# nclass=25)

png("IBD.mantel_corr.ecodist.seplineages.fig5.png",10, 5,units="in",res=300)

palette(brewer.pal(n = 8, name = "Set3"))
cols <- add.alpha(brewer.pal(n = 8, name = "Set2"), 0.8)

#"CA", "WAf","WCA", "DG", "Gab"
par(mar=c(4.1,4.1,.1,.1))
plot(mcc.CA$mgram[,1],mcc.CA$mgram[,3],ylim=c(-.5,+.45),xlab="Distance class (Km)",ylab="",pch=21,col=1,xaxt="n",xlim=c(20,1000),type="n",axes=FALSE,cex=.8,lwd=.7);
axis(1, las=2)
axis(2)
box()
#mtext("  C",side=1,line=-1.7,adj=0,font = 2)
mtext("Mantel correlation",side=2,line=2,cex=0.8)
# ALL
lines(mcc.CA$mgram[,1],mcc.CA$mgram[,3],col=4,lty=1)
for(i in 1:length(mcc.CA$mgram[,1])){
  #if(is.na(mcc.all$mgram[i,4]==TRUE)){print("NA")
  if(mcc.all$mgram[i,2]==0){print("NA")
  }else{points(mcc.CA$mgram[i,1],mcc.CA$mgram[i,3],pch=21,col=4,cex=1,lwd=1)
    arrows(mcc.CA$mgram[i,1],mcc.CA$mgram[i,5],mcc.CA$mgram[i,1], mcc.CA$mgram[i,6], length=0.02, angle=90, code=3,col=4,lty=1,lwd=.6)
    if(mcc.CA$mgram[i,4]<.05){points(mcc.CA$mgram[i,1],mcc.CA$mgram[i,3],pch=21,col=4,bg=4,cex=1,lwd=1)}}}
lines(mcc.WAf$mgram[,1],mcc.WAf$mgram[,3],col=1,lty=1)
for(i in 1:length(mcc.WAf$mgram[,1])){
  #if(is.na(mcc.male$mgram[i,4]==TRUE)){print("NA")
  if(mcc.WAf$mgram[i,2]==0){print("NA")
  }else{points(mcc.WAf$mgram[i,1],mcc.WAf$mgram[i,3],pch=22,col=1,cex=1,lwd=1)
    arrows(mcc.WAf$mgram[i,1],mcc.WAf$mgram[i,5],mcc.WAf$mgram[i,1], mcc.WAf$mgram[i,6], length=0.02, angle=90, code=3,col=1,lty=1,lwd=.6)
    if(mcc.WAf$mgram[i,4]<.05){points(mcc.WAf$mgram[i,1],mcc.WAf$mgram[i,3],pch=22,col=1,bg=1,cex=1,lwd=1)}}}
lines(mcc.WCA$mgram[,1],mcc.WCA$mgram[,3],col=5,lty=1)
for(i in 1:length(mcc.WCA$mgram[,1])){
  #if(is.na(mcc.male$mgram[i,4]==TRUE)){print("NA")
  if(mcc.WCA$mgram[i,2]==0){print("NA")
  }else{points(mcc.WCA$mgram[i,1],mcc.WCA$mgram[i,3],pch=22,col=5,cex=1,lwd=1)
    arrows(mcc.WCA$mgram[i,1],mcc.WCA$mgram[i,5],mcc.WCA$mgram[i,1], mcc.WCA$mgram[i,6], length=0.02, angle=90, code=3,col=5,lty=1,lwd=.6)
    if(mcc.WCA$mgram[i,4]<.05){points(mcc.WCA$mgram[i,1],mcc.WCA$mgram[i,3],pch=22,col=5,bg=5,cex=1,lwd=1)}}}
lines(mcc.DG$mgram[,1],mcc.DG$mgram[,3],col=3,lty=1)
for(i in 1:length(mcc.DG$mgram[,1])){
  #if(is.na(mcc.male$mgram[i,4]==TRUE)){print("NA")
  if(mcc.DG$mgram[i,2]==0){print("NA")
  }else{points(mcc.DG$mgram[i,1],mcc.DG$mgram[i,3],pch=22,col=3,cex=1,lwd=1)
    arrows(mcc.DG$mgram[i,1],mcc.DG$mgram[i,5],mcc.DG$mgram[i,1], mcc.DG$mgram[i,6], length=0.02, angle=90, code=3,col=3,lty=1,lwd=.6)
    if(mcc.DG$mgram[i,4]<.05){points(mcc.DG$mgram[i,1],mcc.DG$mgram[i,3],pch=22,col=3,bg=3,cex=1,lwd=1)}}}
#lines(mcc.Gab$mgram[,1],mcc.Gab$mgram[,3],col=6,lty=1)
abline(h=0,col="red")
legend("topright",legend=c("CA", "WAf","WCA", "DG"),lty=1,pch=1,col=c(4,1,5,3),bty="n")
dev.off()


#######################################################
############ Male & Female Plotting ###################
#######################################################



IBD_dir="/work/idumville/pangolins/RAD/results/IBD/IBD.sex" 
setwd(IBD_dir) #output directory
DGen_orig_auto="/work/idumville/pangolins/RAD/results/IBD/IBD.allsamples/rad.ptri.pango.rad.all_autosomes.allcov_euc_dist.rda"
DGen_orig_xchrom="/work/idumville/pangolins/RAD/results/IBD/IBD.allsamples/rad.ptri.pango.rad.all_xchrom.allcov_euc_dist.rda"
WCA_samples="WCA.txt"
female="newfemales.txt"
male="newmales.txt"
geodata="geodata.txt" #coordinate file of all samples


all_names <- readLines(paste0("../samplelists/", WCA_samples)) #now just R Cameroon
fem_names <- all_names[all_names %in% readLines(paste0("../samplelists/", female))]
male_names <- all_names[all_names %in% readLines(paste0("../samplelists/", male))]


DGen_orig_auto=load(file = DGen_orig_auto)
DGen_orig_auto=cov_euc_dd
DGen_orig_auto <- as.matrix(DGen_orig_auto)
DGen_orig_xchrom=load(file = DGen_orig_xchrom)
DGen_orig_xchrom=cov_euc_dd
DGen_orig_xchrom <- as.matrix(DGen_orig_xchrom)

#Geographic frames
# geographic distances imported as df
Dgeo_orig <- read.table(paste("../samplelists/",geodata,sep=""), row.names = 1)

Dgeo_all <- Dgeo_orig[c(which(rownames(Dgeo_orig) %in% all_names)),]
Dgeo_fem <- Dgeo_orig[c(which(rownames(Dgeo_orig) %in% fem_names)),]
Dgeo_male <- Dgeo_orig[c(which(rownames(Dgeo_orig) %in% male_names)),]

Dgeo_all <- as.dist(spDists(as.matrix(Dgeo_all), longlat=T))
Dgeo_fem <- as.dist(spDists(as.matrix(Dgeo_fem), longlat=T))
Dgeo_male <- as.dist(spDists(as.matrix(Dgeo_male), longlat=T))


#Sorting rownames
DGen_orig_auto_names <- sub("(.*?)(_)", "", x= rownames(DGen_orig_auto))
DGen_orig_auto_names <- str_remove(DGen_orig_auto_names, "Co_") #becasue DR_Co_ kept Co_

DGen_orig_xchrom_names <- sub("(.*?)(_)", "", x= rownames(DGen_orig_xchrom))
DGen_orig_xchrom_names <- str_remove(DGen_orig_xchrom_names, "Co_")#becuase atm concatenated

DGen_all_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% all_names)),c(which(DGen_orig_auto_names %in% all_names))])
DGen_fem_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% fem_names)),c(which(DGen_orig_auto_names %in% fem_names))])
DGen_male_auto <- as.dist(DGen_orig_auto[c(which(DGen_orig_auto_names %in% male_names)),c(which(DGen_orig_auto_names %in% male_names))])

DGen_all_xchrom <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% all_names)),c(which(DGen_orig_xchrom_names %in% all_names))])
DGen_fem_xchrom <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% fem_names)),c(which(DGen_orig_xchrom_names %in% fem_names))])
DGen_male_xchrom <- as.dist(DGen_orig_xchrom[c(which(DGen_orig_xchrom_names %in% male_names)),c(which(DGen_orig_xchrom_names %in% male_names))])


#Estimate the distance at which the IBD is highest ####

source("/work/idumville/pangolins/RAD/bin/R2vsDE_v2.R")
lag=20
minX=50
maxY=0.05

png("IBD.maxdistance.fig1.MF.png",24,10,units="in",res=300)
palette("default")
palette(brewer.pal(n = 8, name = "Set2"))
par(mfrow=c(1,1))
par(mar=c(4.1,4.1,3.1,.1))

#autosomes, all/M/F
mylist<-list(Dgeo=Dgeo_all,Dgen=DGen_all_auto)
Dgen <- DGen_all_auto
Dgeo <- Dgeo_all
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.1,addplot=F,col=1)

mylist<-list(Dgeo=Dgeo_male,Dgen=DGen_male_auto)
Dgen <- DGen_male_auto
Dgeo <- Dgeo_male
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.1,addplot=T,col=2)

mylist<-list(Dgeo=Dgeo_fem,Dgen=DGen_fem_auto)
Dgen <- DGen_fem_auto
Dgeo <- Dgeo_fem
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.1,addplot=T,col=3)

#xchrom all/M/F
mylist<-list(Dgeo=Dgeo_all,Dgen=DGen_all_xchrom)
Dgen <- DGen_all_xchrom
Dgeo <- Dgeo_all
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.1,addplot=T,col=1, lty=2)

mylist<-list(Dgeo=Dgeo_male,Dgen=DGen_male_xchrom)
Dgen <- DGen_male_xchrom
Dgeo <- Dgeo_male
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.1,addplot=T,col=2, lty=2)

mylist<-list(Dgeo=Dgeo_fem,Dgen=DGen_fem_xchrom)
Dgen <- DGen_fem_xchrom
Dgeo <- Dgeo_fem
R2vsDE_1<-R2vsDE(mylist=mylist,fmla=lm(Dgen~Dgeo),lagdist=lag,DEvarname="Dgeo",minX=minX,maxY=0.1,addplot=T,col=3, lty=2)

legend("topleft",lty=c(1,1,1,2,2,2),col=c(1,2,3,1,2,3),legend=c("WCA Autosomes", "WCA M Autosomes", "WCA F Autsomes", "WCA X", "M X", "F X"),title="Gen.Data",cex=0.8)

dev.off()

#Highest distance = largest R2


#Mantel correlograms

#making own classes ; change dgeo and dgen for different populations
mc.all.auto<-mantel.correlog(D.eco=DGen_all_auto,D.geo=Dgeo_all,n.class = 50,nperm = 999)
mc.male.auto<-mantel.correlog(D.eco=DGen_male_auto,D.geo =Dgeo_male,n.class = 50,nperm = 999)
mc.fem.auto<-mantel.correlog(D.eco=DGen_fem_auto,D.geo =Dgeo_fem,n.class = 50,nperm = 999)
mc.all.x<-mantel.correlog(D.eco=DGen_all_xchrom,D.geo =Dgeo_all,n.class = 50,nperm = 999)
mc.male.x<-mantel.correlog(D.eco=DGen_male_xchrom,D.geo =Dgeo_male,n.class = 50,nperm = 999)
mc.fem.x<-mantel.correlog(D.eco=DGen_fem_xchrom,D.geo =Dgeo_fem,n.class = 50,nperm = 999)

png("IBD.mantel_corr.fig2MF.png",20,16,units="in",res=300)
palette("default")
par(mfrow=c(2,3))
par(mar=c(4.1,4.1,.1,.1))
#plot Autosome
plot(mc.all.auto$mantel.res[,1],mc.all.auto$mantel.res[,3],ylim=c(-.2,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,400));
axis(1, las=2)
mtext("  A) Autosome Rural All",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.all.auto$mantel.res[,1],mc.all.auto$mantel.res[,3])
for(i in 1:length(mc.all.auto$mantel.res[,1])){
  if(is.na(mc.all.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.all.auto$mantel.res[i,5]<.05){points(mc.all.auto$mantel.res[i,1],mc.all.auto$mantel.res[i,3],pch=22,col=1,bg=1)}}}
plot(mc.male.auto$mantel.res[,1],mc.male.auto$mantel.res[,3],ylim=c(-.2,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,500));
axis(1, las=2)
mtext("  B) Autosome Males",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.male.auto$mantel.res[,1],mc.male.auto$mantel.res[,3])
for(i in 1:length(mc.male.auto$mantel.res[,1])){
  if(is.na(mc.male.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.male.auto$mantel.res[i,5]<.05){points(mc.male.auto$mantel.res[i,1],mc.male.auto$mantel.res[i,3],pch=22,col=1,bg=1)}}}
plot(mc.fem.auto$mantel.res[,1],mc.fem.auto$mantel.res[,3],ylim=c(-.2,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,500));
axis(1, las=2)
mtext("  C) Autosome Females",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.fem.auto$mantel.res[,1],mc.fem.auto$mantel.res[,3])
for(i in 1:length(mc.fem.auto$mantel.res[,1])){
  if(is.na(mc.fem.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.fem.auto$mantel.res[i,5]<.05){points(mc.fem.auto$mantel.res[i,1],mc.fem.auto$mantel.res[i,3],pch=22,col=1,bg=1)}}}
#plot Xchrom
plot(mc.all.x$mantel.res[,1],mc.all.x$mantel.res[,3],ylim=c(-.2,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,500));
axis(1, las=2)
mtext("  D) Xchrom Rural All",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.all.x$mantel.res[,1],mc.all.x$mantel.res[,3])
for(i in 1:length(mc.all.x$mantel.res[,1])){
  if(is.na(mc.all.x$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.all.x$mantel.res[i,5]<.05){points(mc.all.x$mantel.res[i,1],mc.all.x$mantel.res[i,3],pch=22,col=1,bg=1)}}}
plot(mc.male.x$mantel.res[,1],mc.male.x$mantel.res[,3],ylim=c(-.2,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,500));
axis(1, las=2)
mtext("  E) X Males",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.male.x$mantel.res[,1],mc.male.x$mantel.res[,3])
for(i in 1:length(mc.male.x$mantel.res[,1])){
  if(is.na(mc.male.x$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.male.x$mantel.res[i,5]<.05){points(mc.male.x$mantel.res[i,1],mc.male.x$mantel.res[i,3],pch=22,col=1,bg=1)}}}
plot(mc.fem.x$mantel.res[,1],mc.fem.x$mantel.res[,3],ylim=c(-.2,+.2),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,500));
axis(1, las=2)
mtext("  F) X Females",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.fem.x$mantel.res[,1],mc.fem.x$mantel.res[,3])
for(i in 1:length(mc.fem.x$mantel.res[,1])){
  if(is.na(mc.fem.x$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.fem.x$mantel.res[i,5]<.05){points(mc.fem.x$mantel.res[i,1],mc.fem.auto$mantel.res[i,3],pch=22,col=1,bg=1)}}}
 
dev.off()

#Plot mantel correlograms with M and F and cameroon overlaid

png("IBD.mantel_corr.MF.fig3.png",10,4,units="in",res=300)
palette(brewer.pal(n = 8, name = "Set2"))
par(mfrow=c(1,1))
par(mar=c(4.1,4.1,.1,.1))
#plot(mc.all,ylim=c(-1,1))

#plot F
plot(mc.male.auto$mantel.res[,1],mc.male.auto$mantel.res[,3],ylim=c(-.25,+.25),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,350));
axis(1, las=2)
abline(h=0,col="red")
lines(mc.male.auto$mantel.res[,1],mc.male.auto$mantel.res[,3], col=1, lty=2)
for(i in 1:length(mc.male.auto$mantel.res[,1])){
  if(is.na(mc.male.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.male.auto$mantel.res[i,5]<.05){points(mc.male.auto$mantel.res[i,1],mc.male.auto$mantel.res[i,3],pch=22,col=1,bg=1)}}}

#plot M
points(mc.fem.auto$mantel.res[,1],mc.fem.auto$mantel.res[,3],ylim=c(-.25,+.25),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=23,col=2,xaxt="n",xlim=c(0,300));
axis(1, las=2)
lines(mc.fem.auto$mantel.res[,1],mc.fem.auto$mantel.res[,3], col=2, lty=3)
for(i in 1:length(mc.fem.auto$mantel.res[,1])){
  if(is.na(mc.fem.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.fem.auto$mantel.res[i,5]<.05){points(mc.fem.auto$mantel.res[i,1],mc.fem.auto$mantel.res[i,3],pch=23,col=2,bg=2)}}}  
 
legend("bottomright",legend=c("Male", "Female"),lty=c(2, 3),pch=c(22, 23),col=c(1,2),bty="n")
dev.off()

#######Estimate IBD

ibd.all.auto <- mantel.randtest(DGen_all_auto,Dgeo_all,nrepet = 9999)
ibd.fem.auto <- mantel.randtest(DGen_fem_auto,Dgeo_fem,nrepet = 9999)
ibd.male.auto <- mantel.randtest(DGen_male_auto,Dgeo_male,nrepet = 9999)

###overwriting BE CAREFUL
ibd.all.auto <- mantel.randtest(DGen_all_xchrom,Dgeo_all,nrepet = 9999)
ibd.fem.auto <- mantel.randtest(DGen_fem_xchrom,Dgeo_fem,nrepet = 9999)
ibd.male.auto <- mantel.randtest(DGen_male_xchrom,Dgeo_male,nrepet = 9999)


ibd.fem.x <- mantel.randtest(DGen_fem_xchrom,Dgeo_fem,nrepet = 9999)

#### Final figure IBD
#plotting autosome urban, rural and cameroon
png("IBD.main_fig_mantel_corr.fig4_MF_Xchrom.png",10,5,units="in",res=300)
palette("default")
palette(brewer.pal(n = 8, name = "Set2"))
cols <- add.alpha(brewer.pal(n = 8, name = "Set2"), 0.5)

par(mfrow=c(1,1))
par(mar=c(4.1,5.1,.1,.1))

# IBD corr + hist
X1<-0.7*max(Dgeo_all)
X<-0.8*max(Dgeo_all)
X2<-0.6*max(Dgeo_all)
X3<-0.9*max(Dgeo_all)
plot(Dgeo_all, DGen_all_auto, pch=21,cex=.6, col="grey20", bg=cols[1], xlab ="Geographic distance (Km)", lwd=0.2, #cols[3]
    ylab="Genetic distance",main="",xaxt="n")
axis(side = 1,las=3,cex=0.9)
abline(lm(DGen_all_auto~Dgeo_all), lty=2,col=1) #cols[3]

points(Dgeo_male, DGen_male_auto, pch=21, cex=.6, col="grey20", bg=2, lwd=0.2)
abline(lm(DGen_male_auto~Dgeo_male), lty=3,col=2)

points(Dgeo_fem, DGen_fem_auto, pch=21, cex=.6, col="grey20", bg=3, lwd=0.2)
abline(lm(DGen_fem_auto~Dgeo_fem), lty=3,col=3)

text(X2,0.225*max(DGen_all_auto),"All Auto", cex=.6)
text(X2,0.175*max(DGen_all_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.all.auto$obs)^2,3))),cex=.6,col="black")
text(X2,0.125*max(DGen_all_auto),substitute(italic(p)<=rr,list(rr=round(ibd.all.auto$pvalue,4))),cex=.6,col="black")

text(X1,0.225*max(DGen_all_auto),"Male Auto", cex=.6)
text(X1,0.175*max(DGen_all_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.male.auto$obs)^2,3))),cex=.6,col="black")
text(X1,0.125*max(DGen_all_auto),substitute(italic(p)<=rr,list(rr=round(ibd.male.auto$pvalue,4))),cex=.6,col="black")

text(X,0.225*max(DGen_all_auto),"Female Auto", cex=.6)
text(X,0.175*max(DGen_all_auto),substitute(italic(R)^2==rr,list(rr=round((ibd.fem.auto$obs)^2,3))),cex=.6,col="black")
text(X,0.125*max(DGen_all_auto),substitute(italic(p)<=rr,list(rr=round(ibd.fem.auto$pvalue,4))),cex=.6,col="black")

legend("topright",pch=21,col="grey20", pt.bg=c(1,2,3,4),legend=c("WCA Autosomes", "Male Autosomes", "Female Autosomes"),cex=0.8)

dev.off()


## Ecodist Figure

library(ecodist)

#also done with xchrom


for (val in c(30,50,65,75))
{
stepsize=val
mcc.all<-mgram(species.d=DGen_all_auto, space.d=Dgeo_all, stepsize=stepsize)# nclass=25)
mcc.female<-mgram(species.d=DGen_fem_auto, space.d=Dgeo_fem, stepsize=stepsize)# nclass=25)
mcc.male<-mgram(species.d=DGen_male_auto, space.d=Dgeo_male, stepsize=stepsize)# nclass=25)
mcc.female.x<-mgram(species.d=DGen_fem_xchrom, space.d=Dgeo_fem, stepsize=stepsize)


png(paste0("IBD.mantel_corr.ecodist.MF_Autosomes_andX_step", stepsize, ".fig5.png"),10,4,units="in",res=300)


palette(brewer.pal(n = 8, name = "Set2"))
cols <- add.alpha(brewer.pal(n = 8, name = "Set2"), 0.5)

par(mar=c(4.1,4.1,.1,0.1))
plot(mcc.all$mgram[,1],mcc.all$mgram[,3],ylim=c(-.25,+.25),xlab="Distance class (Km)",ylab="",pch=21,col=1,xaxt="n",xlim=c(20,700),type="n",axes=FALSE,cex=.8,lwd=.7);
axis(1, las=2)
axis(2)
box()
#mtext("  C",side=1,line=-1.7,adj=0,font = 2)
mtext("Mantel correlation",side=2,line=2,cex=0.8)
# ALL
lines(mcc.all$mgram[,1],mcc.all$mgram[,3],col=1,lty=1)
for(i in 1:length(mcc.all$mgram[,1])){
  #if(is.na(mcc.all$mgram[i,4]==TRUE)){print("NA")
  if(mcc.all$mgram[i,2]==0){print("NA")
  }else{points(mcc.all$mgram[i,1],mcc.all$mgram[i,3],pch=21,col=1,cex=1,lwd=1)
    arrows(mcc.all$mgram[i,1],mcc.all$mgram[i,5],mcc.all$mgram[i,1], mcc.all$mgram[i,6], length=0.02, angle=90, code=3,col=1,lty=1,lwd=.6)
    if(mcc.all$mgram[i,4]<.05){points(mcc.all$mgram[i,1],mcc.all$mgram[i,3],pch=21,col=1,bg=1,cex=1,lwd=1)}}}
# MALES
lines(mcc.male$mgram[,1],mcc.male$mgram[,3],col=2,lty=2)
for(i in 1:length(mcc.male$mgram[,1])){
  #if(is.na(mcc.male$mgram[i,4]==TRUE)){print("NA")
  if(mcc.male$mgram[i,2]==0){print("NA")
  }else{points(mcc.male$mgram[i,1],mcc.male$mgram[i,3],pch=22,col=2,cex=1,lwd=1)
    arrows(mcc.male$mgram[i,1],mcc.male$mgram[i,5],mcc.male$mgram[i,1], mcc.male$mgram[i,6], length=0.02, angle=90, code=3,col=2,lty=1,lwd=.6)
    if(mcc.male$mgram[i,4]<.05){points(mcc.male$mgram[i,1],mcc.male$mgram[i,3],pch=22,col=2,bg=2,cex=1,lwd=1)}}}
# Females
lines(mcc.female$mgram[,1],mcc.female$mgram[,3],col=3,lty=3)
for(i in 1:length(mcc.female$mgram[,1])){
  #if(is.na(mcc.male$mgram[i,4]==TRUE)){print("NA")
  if(mcc.female$mgram[i,2]==0){print("NA")
  }else{points(mcc.female$mgram[i,1],mcc.female$mgram[i,3],pch=23,col=3,cex=1,lwd=1)
        arrows(mcc.female$mgram[i,1],mcc.female$mgram[i,5],mcc.female$mgram[i,1], mcc.female$mgram[i,6], length=0.03, angle=90, code=3,col=3,lty=1,lwd=.6)
    if(mcc.female$mgram[i,4]<.05){points(mcc.female$mgram[i,1],mcc.female$mgram[i,3],pch=23,col=3,bg=3,cex=1,lwd=1)}}}
# Females X
lines(mcc.female.x$mgram[,1],mcc.female.x$mgram[,3],col=4,lty=4)
for(i in 1:length(mcc.female.x$mgram[,1])){
  #if(is.na(mcc.male$mgram[i,4]==TRUE)){print("NA")
  if(mcc.female.x$mgram[i,2]==0){print("NA")
  }else{points(mcc.female.x$mgram[i,1],mcc.female.x$mgram[i,3],pch=24,col=4,cex=1,lwd=1)
        arrows(mcc.female.x$mgram[i,1],mcc.female.x$mgram[i,5],mcc.female.x$mgram[i,1], mcc.female.x$mgram[i,6], length=0.03, angle=90, code=3,col=4,lty=1,lwd=.6)
    if(mcc.female.x$mgram[i,4]<.05){points(mcc.female.x$mgram[i,1],mcc.female.x$mgram[i,3],pch=24,col=4,bg=4,cex=1,lwd=1)}}}
abline(h=0,col="red")
legend("topright",legend=c("All Autosomes","Male Autosomes","Female Autosomes", "Female X"),lty=c(1,2,3,4),pch=c(21,22,23,24),col=c(1,2,3, 4),bty="n")
dev.off()
}


####################   END MALE AND FEMALE PLOTTING




###### Extracting values of matrix above threshold
Dgen_concat_cameroon
Dgeo_cameroon

#mat1[which(mat2 == 0)] <- NA

fit <- lm(DGen_concat_cameroon ~ Dgeo_cameroon)
fit
message("manually noted in script")
#doing manually

temp <- as.matrix(Dgeo_cameroon)
modelled_DGen_concat_cameroon <- apply(temp, c(1,2), FUN = function(x){(0.0004182*x) + 0.6156028}) #add 0.15 to intercept
modelled_DGen_concat_cameroon <- as.dist(modelled_DGen_concat_cameroon)

equatedmatix <- apply(geodistance, (lm(DGen_concat~Dgeo) + x))
DGen_concat_cameroon <- as.matrix(DGen_concat_cameroon)
modelled_DGen_concat_cameroon <- as.matrix(modelled_DGen_concat_cameroon)
which(DGen_concat_cameroon > modelled_DGen_concat_cameroon)
newmatrix <- as.dist(DGen_concat_cameroon[which(DGen_concat_cameroon > modelled_DGen_concat_cameroon)])

DGen_concat_cameroon[which(DGen_concat_cameroon < modelled_DGen_concat_cameroon)] <- NA
DGen_concat_cameroon <- as.dist(DGen_concat_cameroon)

##Replotting the kept values to check kept the right ones
png("IBD.main_fig_mantel_corr.concatcameroon.ed.fig.png",16,8,units="in",res=300)
palette(brewer.pal(n = 8, name = "Set3"))
cols <- add.alpha(brewer.pal(n = 8, name = "Set3"), 0.5)
par(mfrow=c(1,1))
par(mar=c(4.1,5.1,.1,.1))
# IBD corr + hist
X<-0.9*max(Dgeo)
X2<-0.65*max(Dgeo)
plot(Dgeo, DGen_concat, pch=18,cex=.6, col=cols[3], xlab ="Geographic distance (Km)",
     ylab="Genetic distance",main="",xaxt="n")
axis(side = 1,las=3,cex=0.9)
abline(lm(DGen_concat~Dgeo), lty=2,col=3)
points(Dgeo_cameroon, DGen_concat_cameroon, pch=18, cex=.6, col=cols[4])
abline(lm(DGen_concat_cameroon~Dgeo_cameroon), lty=3,col=4)
points(Dgeo_cameroon, modelled_DGen_concat_cameroon, pch=18, cex=.6, col=cols[5])
abline(lm(modelled_DGen_concat_cameroon~Dgeo_cameroon), lty=3,col=5)
text(X,0.1*max(DGen_concat),"All", cex=.8)
text(X,0.0,substitute(italic(R)^2==rr,list(rr=round((ibd.concat$obs)^2,3))),cex=.8,col="black")
text(X,0.05*max(DGen_concat),substitute(italic(p)<=rr,list(rr=round(ibd.concat$pvalue,4))),cex=.8,col="black")
text(X2,0.1*max(DGen_concat),"Cameroon", cex=.8)
text(X2,0.0,substitute(italic(R)^2==rr,list(rr=round((ibd.cameroon$obs)^2,3))),cex=.8,col="black")
text(X2,0.05*max(DGen_concat),substitute(italic(p)<=rr,list(rr=round(ibd.cameroon$pvalue,4))),cex=.8,col="black")
legend("topright",pch=c(18,18,18),col=c(3,4,5),legend=c("All Concatenated Samples","Concatenated Cameroon", "Modelled Cameroon"),title="Dataset",cex=0.8)
dev.off()


#now removing any rows and columsn where just NA
DGen_concat_cameroon <- as.matrix(DGen_concat_cameroon)

DGen_concat_cameroon <-DGen_concat_cameroon[rowSums(is.na(DGen_concat_cameroon )) != ncol(DGen_concat_cameroon ), ]
DGen_concat_cameroon <-DGen_concat_cameroon[colSums(is.na(DGen_concat_cameroon )) != nrow(DGen_concat_cameroon ), ]

write.table(DGen_concat_cameroon, file = "./Cameroon_extracted_biggerGD.txt", sep="\t")


q(save="no")







##########SCRAP



#plot modern hist ecodist mantel correlogram
par(mar=c(.1,.1,.1,4.1))
plot(mcc.all$mgram[,1],mcc.all$mgram[,3],ylim=c(-.46,+.3),xlab="Distance class index (in Km)",ylab=,pch=21,col=1,xaxt="n",xlim=c(35,5900),type="n",
     axes=FALSE,cex=0.6);
axis(4)
box()
mtext("  B",side=1,line=-1.7,adj=0,font = 2)
mtext("Mantel correlation",side=4,line=3,cex=0.6)
# ALL
lines(mcc.all$mgram[,1],mcc.all$mgram[,3],col=1,lty=1)
for(i in 1:length(mcc.all$mgram[,1])){
  if(mcc.all$mgram[i,2]==0){print("NA")
  }else{points(mcc.all$mgram[i,1],mcc.all$mgram[i,3],pch=21,col=1,cex=.8,lwd=.7)
    arrows(mcc.all$mgram[i,1],mcc.all$mgram[i,5],mcc.all$mgram[i,1], mcc.all$mgram[i,6], length=0.02, angle=90, code=3,col=1,lty=1,lwd=.6)
    if(mcc.all$mgram[i,4]<.05){points(mcc.all$mgram[i,1],mcc.all$mgram[i,3],pch=21,col=1,bg=1,cex=.8,lwd=.7)}}}
# Modern
lines(mcc.m$mgram[,1],mcc.m$mgram[,3],col=2,lty=2)
for(i in 1:length(mcc.m$mgram[,1])){
  if(mcc.m$mgram[i,2]==0){print("NA")
  }else{points(mcc.m$mgram[i,1]-20,mcc.m$mgram[i,3],pch=22,col=2,cex=.8,lwd=.7)
    arrows(mcc.m$mgram[i,1]-20,mcc.m$mgram[i,5],mcc.m$mgram[i,1]-20, mcc.m$mgram[i,6], length=0.02, angle=90, code=3,col=2,lty=1,lwd=.6)
    if(mcc.m$mgram[i,4]<.05){points(mcc.m$mgram[i,1]-20,mcc.m$mgram[i,3],pch=22,col=2,bg=2,cex=.8,lwd=.7)}}}
# Historic
lines(mcc.h$mgram[,1]+20,mcc.h$mgram[,3],col=3,lty=3)
for(i in 1:length(mcc.h$mgram[,1])){
  if(mcc.h$mgram[i,2]==0){print("NA")
  }else{points(mcc.h$mgram[i,1]+20,mcc.h$mgram[i,3],pch=23,col=3,cex=.8,lwd=.7)
        arrows(mcc.h$mgram[i,1]+20,mcc.h$mgram[i,5],mcc.h$mgram[i,1]+20, mcc.h$mgram[i,6], length=0.02, angle=90, code=3,col=3,lty=1,lwd=.6)
    if(mcc.h$mgram[i,4]<.05){points(mcc.h$mgram[i,1]+20,mcc.h$mgram[i,3],pch=23,col=3,bg=3,cex=.8,lwd=.7)}}}
abline(h=0,col="red")
legend("topright",legend=c("All","Modern","Historic"),lty=c(1,2,3),pch=c(21,22,23),col=c(1,2,3),bty="n")

# BY SEX
par(mar=c(4.1,.1,.1,4.1))
plot(mcc.all$mgram[,1],mcc.all$mgram[,3],ylim=c(-.46,+.3),xlab="Distance class (Km)",ylab="Mantel correlation",pch=21,col=1,xaxt="n",xlim=c(35,5900),type="n",axes=FALSE,cex=.8,lwd=.7);
axis(1, las=2)
axis(4)
box()
mtext("  C",side=1,line=-1.7,adj=0,font = 2)
mtext("Mantel correlation",side=4,line=3,cex=0.6)
# ALL
lines(mcc.all$mgram[,1],mcc.all$mgram[,3],col=1,lty=1)
for(i in 1:length(mcc.all$mgram[,1])){
  #if(is.na(mcc.all$mgram[i,4]==TRUE)){print("NA")
  if(mcc.all$mgram[i,2]==0){print("NA")
  }else{points(mcc.all$mgram[i,1],mcc.all$mgram[i,3],pch=21,col=1,cex=.8,lwd=.7)
    arrows(mcc.all$mgram[i,1],mcc.all$mgram[i,5],mcc.all$mgram[i,1], mcc.all$mgram[i,6], length=0.02, angle=90, code=3,col=1,lty=1,lwd=.6)
    if(mcc.all$mgram[i,4]<.05){points(mcc.all$mgram[i,1],mcc.all$mgram[i,3],pch=21,col=1,bg=1,cex=.8,lwd=.7)}}}
# MALES
lines(mcc.male$mgram[,1],mcc.male$mgram[,3],col=2,lty=2)
for(i in 1:length(mcc.male$mgram[,1])){
  #if(is.na(mcc.male$mgram[i,4]==TRUE)){print("NA")
  if(mcc.male$mgram[i,2]==0){print("NA")
  }else{points(mcc.male$mgram[i,1]-20,mcc.male$mgram[i,3],pch=22,col=2,cex=.8,lwd=.7)
    arrows(mcc.male$mgram[i,1]-20,mcc.male$mgram[i,5],mcc.male$mgram[i,1]-20, mcc.male$mgram[i,6], length=0.02, angle=90, code=3,col=2,lty=1,lwd=.6)
    if(mcc.male$mgram[i,4]<.05){points(mcc.male$mgram[i,1]-20,mcc.male$mgram[i,3],pch=22,col=2,bg=2,cex=.8,lwd=.7)}}}
# Females
lines(mcc.female$mgram[,1]+20,mcc.female$mgram[,3],col=3,lty=3)
for(i in 1:length(mcc.female$mgram[,1])){
  #if(is.na(mcc.male$mgram[i,4]==TRUE)){print("NA")
  if(mcc.female$mgram[i,2]==0){print("NA")
  }else{points(mcc.female$mgram[i,1]+20,mcc.female$mgram[i,3],pch=23,col=3,cex=.8,lwd=.7)
        arrows(mcc.female$mgram[i,1]+20,mcc.female$mgram[i,5],mcc.female$mgram[i,1]+20, mcc.female$mgram[i,6], length=0.03, angle=90, code=3,col=3,lty=1,lwd=.6)
    if(mcc.female$mgram[i,4]<.05){points(mcc.female$mgram[i,1]+20,mcc.female$mgram[i,3],pch=23,col=3,bg=3,cex=.8,lwd=.7)}}}
abline(h=0,col="red")
legend("topright",legend=c("All","Males","Females"),lty=c(1,2,3),pch=c(21,22,23),col=c(1,2,3),bty="n")
dev.off()


png("IBD.mantel_corr.AMF.fig.png",10,12,units="in",res=300)

par(mfrow=c(3,1))
par(mar=c(4.1,4.1,.1,.1))
#plot(mc.all,ylim=c(-1,1))
plot(mc.all$mantel.res[,1],mc.all$mantel.res[,3],ylim=c(-.25,+.65),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,3000));
axis(1, las=2)
mtext("  A) All samples",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.all$mantel.res[,1],mc.all$mantel.res[,3])
for(i in 1:length(mc.all$mantel.res[,1])){
  if(is.na(mc.all$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.all$mantel.res[i,5]<.05){points(mc.all$mantel.res[i,1],mc.all$mantel.res[i,3],pch=22,col=1,bg=1)}}}
plot(mc.male$mantel.res[,1],mc.male$mantel.res[,3],ylim=c(-.25,+.65),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,3000));
axis(1, las=2)
mtext("  B) Male samples",side=3,line=-1.7,adj=0,font = 2)
abline(h=0,col="red")
lines(mc.male$mantel.res[,1],mc.male$mantel.res[,3])
for(i in 1:length(mc.male$mantel.res[,1])){
  if(is.na(mc.male$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.male$mantel.res[i,5]<.05){points(mc.male$mantel.res[i,1],mc.male$mantel.res[i,3],pch=22,col=1,bg=1)}}}
# plot historic
plot(mc.female$mantel.res[,1],mc.female$mantel.res[,3],ylim=c(-.25,+.65),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=1,xaxt="n",xlim=c(0,3000));
mtext("  C) Female samples",side=3,line=-1.7,adj=0,font = 2)
axis(1, las=2)
abline(h=0,col="red")
lines(mc.female$mantel.res[,1],mc.female$mantel.res[,3])
for(i in 1:length(mc.female$mantel.res[,1])){
  if(is.na(mc.female$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.female$mantel.res[i,5]<.05){points(mc.female$mantel.res[i,1],mc.female$mantel.res[i,3],pch=22,col=1,bg=1)}}}

dev.off()
```
Plot mantel correlograms males, females





```
Mantel correlogram with ecodist male and female
```{r}
library(ecodist)
stepsize=100
mcc.a<-mgram(species.d=DGen_concat, space.d=Dgeo, stepsize=stepsize)# nclass=25)
mcc.cam<-mgram(species.d=DGen_concat_cameroon, space.d=Dgeo_cameroon, stepsize=stepsize)# nclass=25)


```
plot modern hist ecodist mantel correlogram
```{r}
png("IBD.mantel_corr.AMH.AMF.ed.fig.png",10,8,units="in",res=300)

par(mfrow=c(2,1))
par(mar=c(4.1,4.1,.1,.1))
plot(mcc.all$mgram[,1],mcc.all$mgram[,3],ylim=c(-.45,+.45),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=21,col=1,xaxt="n",xlim=c(35,5900),type="n");
axis(1, las=2)
mtext("  A) Temporal analysis",side=3,line=-1.7,adj=0,font = 2)

lines(mcc.all$mgram[,1],mcc.all$mgram[,3],col=1,lty=1)
for(i in 1:length(mcc.all$mgram[,1])){
  #if(is.na(mcc.all$mgram[i,4]==TRUE)){print("NA")
  if(mcc.all$mgram[i,2]==0){print("NA")
  }else{points(mcc.all$mgram[i,1],mcc.all$mgram[i,3],pch=21,col=1)
    arrows(mcc.all$mgram[i,1],mcc.all$mgram[i,5],mcc.all$mgram[i,1], mcc.all$mgram[i,6], length=0.04, angle=90, code=3,col=1,lty=1)
    if(mcc.all$mgram[i,4]<.05){points(mcc.all$mgram[i,1],mcc.all$mgram[i,3],pch=21,col=1,bg=1)}}}

#plot(mc.all,ylim=c(-1,1))
#plot(mcc.male$mgram[,1],mcc.male$mgram[,3],ylim=c(-.45,+.45),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=2,xaxt="n",xlim=c(0,5200),type="n");
#axis(1, las=2)
#mtext("  B) Male samples",side=3,line=-1.7,adj=0,font = 2)

lines(mcc.m$mgram[,1],mcc.m$mgram[,3],col=2,lty=2)
for(i in 1:length(mcc.m$mgram[,1])){
  #if(is.na(mcc.m$mgram[i,4]==TRUE)){print("NA")
  if(mcc.m$mgram[i,2]==0){print("NA")
  }else{points(mcc.m$mgram[i,1]-20,mcc.m$mgram[i,3],pch=22,col=2)
    arrows(mcc.m$mgram[i,1]-20,mcc.m$mgram[i,5],mcc.m$mgram[i,1]-20, mcc.m$mgram[i,6], length=0.04, angle=90, code=3,col=2,lty=1)
    if(mcc.m$mgram[i,4]<.05){points(mcc.m$mgram[i,1]-20,mcc.m$mgram[i,3],pch=22,col=2,bg=2)}}}

# plot females
#points(mcc.male$mgram[,1],mcc.male$mgram[,3],ylim=c(-.42,+.42),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=2,xaxt="n",xlim=c(0,5200),type="n");
#axis(1, las=2)
#mtext("  B) Male samples",side=3,line=-1.7,adj=0,font = 2)
#abline(h=0,col="red")
lines(mcc.h$mgram[,1]+20,mcc.h$mgram[,3],col=3,lty=3)
for(i in 1:length(mcc.h$mgram[,1])){
  #if(is.na(mcc.male$mgram[i,4]==TRUE)){print("NA")
  if(mcc.h$mgram[i,2]==0){print("NA")
  }else{points(mcc.h$mgram[i,1]+20,mcc.h$mgram[i,3],pch=23,col=3)
        arrows(mcc.h$mgram[i,1]+20,mcc.h$mgram[i,5],mcc.h$mgram[i,1]+20, mcc.h$mgram[i,6], length=0.04, angle=90, code=3,col=3,lty=1)
    if(mcc.h$mgram[i,4]<.05){points(mcc.h$mgram[i,1]+20,mcc.h$mgram[i,3],pch=23,col=3,bg=3)}}}
abline(h=0,col="red")
legend("topright",legend=c("All","Modern","Historic"),lty=c(1,2,3),pch=c(21,22,23),col=c(1,2,3),bty="n")
#dev.off()

#png("IBD.mantel_corr.AMF.ed.fig.png",10,4,units="in",res=300)

#par(mfrow=c(1,1))
#par(mar=c(4.1,4.1,.1,.1))
plot(mcc.all$mgram[,1],mcc.all$mgram[,3],ylim=c(-.45,+.45),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=21,col=1,xaxt="n",xlim=c(35,5900),type="n");
axis(1, las=2)
mtext("  B) Analysis by sex",side=3,line=-1.7,adj=0,font = 2)

lines(mcc.all$mgram[,1],mcc.all$mgram[,3],col=1,lty=1)
for(i in 1:length(mcc.all$mgram[,1])){
  #if(is.na(mcc.all$mgram[i,4]==TRUE)){print("NA")
  if(mcc.all$mgram[i,2]==0){print("NA")
  }else{points(mcc.all$mgram[i,1],mcc.all$mgram[i,3],pch=21,col=1)
    arrows(mcc.all$mgram[i,1],mcc.all$mgram[i,5],mcc.all$mgram[i,1], mcc.all$mgram[i,6], length=0.04, angle=90, code=3,col=1,lty=1)
    if(mcc.all$mgram[i,4]<.05){points(mcc.all$mgram[i,1],mcc.all$mgram[i,3],pch=21,col=1,bg=1)}}}

#plot(mc.all,ylim=c(-1,1))
#plot(mcc.male$mgram[,1],mcc.male$mgram[,3],ylim=c(-.45,+.45),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=2,xaxt="n",xlim=c(0,5200),type="n");
#axis(1, las=2)
#mtext("  B) Male samples",side=3,line=-1.7,adj=0,font = 2)

lines(mcc.male$mgram[,1],mcc.male$mgram[,3],col=2,lty=2)
for(i in 1:length(mcc.male$mgram[,1])){
  #if(is.na(mcc.male$mgram[i,4]==TRUE)){print("NA")
  if(mcc.male$mgram[i,2]==0){print("NA")
  }else{points(mcc.male$mgram[i,1]-20,mcc.male$mgram[i,3],pch=22,col=2)
    arrows(mcc.male$mgram[i,1]-20,mcc.male$mgram[i,5],mcc.male$mgram[i,1]-20, mcc.male$mgram[i,6], length=0.04, angle=90, code=3,col=2,lty=1)
    if(mcc.male$mgram[i,4]<.05){points(mcc.male$mgram[i,1]-20,mcc.male$mgram[i,3],pch=22,col=2,bg=2)}}}
# plot females
#points(mcc.male$mgram[,1],mcc.male$mgram[,3],ylim=c(-.42,+.42),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=2,xaxt="n",xlim=c(0,5200),type="n");
#axis(1, las=2)
#mtext("  B) Male samples",side=3,line=-1.7,adj=0,font = 2)
#abline(h=0,col="red")
lines(mcc.female$mgram[,1]+20,mcc.female$mgram[,3],col=3,lty=3)
for(i in 1:length(mcc.female$mgram[,1])){
  #if(is.na(mcc.male$mgram[i,4]==TRUE)){print("NA")
  if(mcc.female$mgram[i,2]==0){print("NA")
  }else{points(mcc.female$mgram[i,1]+20,mcc.female$mgram[i,3],pch=23,col=3)
        arrows(mcc.female$mgram[i,1]+20,mcc.female$mgram[i,5],mcc.female$mgram[i,1]+20, mcc.female$mgram[i,6], length=0.04, angle=90, code=3,col=3,lty=1)
    if(mcc.female$mgram[i,4]<.05){points(mcc.female$mgram[i,1]+20,mcc.female$mgram[i,3],pch=23,col=3,bg=3)}}}
abline(h=0,col="red")
legend("topright",legend=c("All","Males","Females"),lty=c(1,2,3),pch=c(21,22,23),col=c(1,2,3),bty="n")
dev.off()
```

plot ecodist mantel correlogram adults
```{r}
png("IBD.mantel_corr.MF.ad.ed.fig.png",10,4,units="in",res=300)

par(mfrow=c(1,1))
par(mar=c(4.1,4.1,.1,.1))
#plot(mc.all,ylim=c(-1,1))
plot(mcc.male_ad$mgram[,1],mcc.male_ad$mgram[,3],ylim=c(-.45,+.45),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=2,xaxt="n",xlim=c(0,5200),type="n");
axis(1, las=2)
#mtext("  B) male_ad samples",side=3,line=-1.7,adj=0,font = 2)

lines(mcc.male_ad$mgram[,1],mcc.male_ad$mgram[,3],col=2,lty=2)
for(i in 1:length(mcc.male_ad$mgram[,1])){
  #if(is.na(mcc.male_ad$mgram[i,4]==TRUE)){print("NA")
  if(mcc.male_ad$mgram[i,2]==0){print("NA")
  }else{points(mcc.male_ad$mgram[i,1],mcc.male_ad$mgram[i,3],pch=22,col=2)
    arrows(mcc.male_ad$mgram[i,1],mcc.male_ad$mgram[i,5],mcc.male_ad$mgram[i,1], mcc.male_ad$mgram[i,6], length=0.04, angle=90, code=3,col=2,lty=1)
    if(mcc.male_ad$mgram[i,4]<.05){points(mcc.male_ad$mgram[i,1],mcc.male_ad$mgram[i,3],pch=22,col=2,bg=2)}}}
# plot female_ads
#points(mcc.male_ad$mgram[,1],mcc.male_ad$mgram[,3],ylim=c(-.42,+.42),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=2,xaxt="n",xlim=c(0,5200),type="n");
#axis(1, las=2)
#mtext("  B) male_ad samples",side=3,line=-1.7,adj=0,font = 2)
#abline(h=0,col="red")
lines(mcc.female_ad$mgram[,1],mcc.female_ad$mgram[,3],col=3,lty=3)
for(i in 1:length(mcc.female_ad$mgram[,1])){
  #if(is.na(mcc.male_ad$mgram[i,4]==TRUE)){print("NA")
  if(mcc.female_ad$mgram[i,2]==0){print("NA")
  }else{points(mcc.female_ad$mgram[i,1],mcc.female_ad$mgram[i,3],pch=23,col=3)
        arrows(mcc.female_ad$mgram[i,1],mcc.female_ad$mgram[i,5],mcc.female_ad$mgram[i,1], mcc.female_ad$mgram[i,6], length=0.04, angle=90, code=3,col=3,lty=1)
    if(mcc.female_ad$mgram[i,4]<.05){points(mcc.female_ad$mgram[i,1],mcc.female_ad$mgram[i,3],pch=23,col=3,bg=3)}}}
abline(h=0,col="red")
legend("topright",legend=c("Males","Females"),lty=c(2,3),pch=c(22,23),col=c(2,3),bty="n")
dev.off()
```



#Old scrap
loc="autosomes" #xchrom #autosomes
prefix="lowreads"
fileprefix="ruralhighread_wobio" #gab

DGen_concat="/work/idumville/pangolins/RAD/results/pcangsd.concat.ruralhighread_wobio/ruralhighread_wobio.pango.rad.autosomes.allcov_euc_dist.rda"
DGen_x="/work/idumville/pangolins/RAD/results/pcangsd.xchrom.ruralhighread_wobio/ruralhighread_wobio.pango.rad.xchrom.allcov_euc_dist.rda"
DGen_auto="/work/idumville/pangolins/RAD/results/pcangsd.autosomes.ruralhighread_wobio/ruralhighread_wobio.pango.rad.autosomes.allcov_euc_dist.rda"
Cam_samples="/work/idumville/pangolins/RAD/results/IBD.ruralhighread_wobio/cameroon.txt"
samplelist="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads/autosomesruralhighread_wobiosamples.txt"
geodata="geodataruralhighread_wobio.txt"

#gendist="rad.ptri.pango.rad.all_xchrom.allcov_euc_dist.rda"
#running Rscript
module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; Rscript $bin/IBD_plotting.R $IBD_dir geodata${fileprefix}.txt $samplelist $DGen_concat $DGen_x $DGen_auto $Cam_samples



samplelist=readLines(samplelist)
DGen_concat=load(file = DGen_concat)
DGen_concat=cov_euc_dd #all R objects saved with same name
DGen_x=load(file = DGen_x)
DGen_x=cov_euc_dd
DGen_auto=load(file = DGen_auto)
DGen_auto=cov_euc_dd

#making cameron matrixes; list of which samples in sample list = cameroon
camer <- readLines(Cam_samples)
cameroonnums <- which(samplelist %in% camer)
cameroonnames <- camer[which(camer %in% samplelist)]






#changing row and col names, creating cameroon objects (change to matrix for ease of subsetting)
Dgeo <- as.matrix(Dgeo)
camDgeo <- Dgeo[c(cameroonnums), c(cameroonnums)]
colnames(camDgeo) <- cameroonnames
rownames(camDgeo) <- cameroonnames
DGen_concat <- as.matrix(DGen_concat)
DGen_concat_cameroon <- DGen_concat[c(cameroonnums), c(cameroonnums)]
colnames(Dgeo) <- rownames(DGen_concat)
rownames(Dgeo) <- rownames(DGen_concat)
Dgeo <- as.dist(Dgeo)
Dgeo_cameroon <- as.dist(camDgeo)
DGen_concat <- as.dist(DGen_concat)
DGen_concat_cameroon <- as.dist(DGen_concat_cameroon)
##making cameroon genetic distance



#plot cam autosomes
points(mc.cam.auto$mantel.res[,1],mc.cam.auto$mantel.res[,3],ylim=c(-.25,+.65),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=3,xaxt="n",xlim=c(0,3000));
axis(1, las=2)
lines(mc.cam.auto$mantel.res[,1],mc.cam.auto$mantel.res[,3], col=4, lty=3)
for(i in 1:length(mc.cam.auto$mantel.res[,1])){
  if(is.na(mc.cam.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.cam.auto$mantel.res[i,5]<.05){points(mc.cam.auto$mantel.res[i,1],mc.cam.auto$mantel.res[i,3],pch=22,col=4,bg=4)}}}
dec
#plot all autosomes
plot(mc.all.auto$mantel.res[,1],mc.all.auto$mantel.res[,3],ylim=c(-.2,+.15),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=2,xaxt="n",xlim=c(0,400));
axis(1, las=2)
abline(h=0,col="red")
lines(mc.all.auto$mantel.res[,1],mc.all.auto$mantel.res[,3], col=2, lty=1)
for(i in 1:length(mc.all.auto$mantel.res[,1])){
  if(is.na(mc.all.auto$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.all.auto$mantel.res[i,5]<.05){points(mc.all.auto$mantel.res[i,1],mc.all.auto$mantel.res[i,3],pch=22,col=2,bg=2)}}}
#plot all x
points(mc.all.x$mantel.res[,1],mc.all.x$mantel.res[,3],ylim=c(-.25,+.65),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=3,xaxt="n",xlim=c(0,3000));
axis(1, las=2)
lines(mc.all.x$mantel.res[,1],mc.all.x$mantel.res[,3], col=5, lty=1)
for(i in 1:length(mc.all.x$mantel.res[,1])){
  if(is.na(mc.all.x$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.all.x$mantel.res[i,5]<.05){points(mc.all.x$mantel.res[i,1],mc.all.x$mantel.res[i,3],pch=22,col=5,bg=5)}}} 
#plot cam x
points(mc.cam.x$mantel.res[,1],mc.cam.x$mantel.res[,3],ylim=c(-.25,+.65),xlab="Distance class index (in Km.)",ylab="Mantel correlation",pch=22,col=3,xaxt="n",xlim=c(0,3000));
axis(1, las=2)
lines(mc.cam.x$mantel.res[,1],mc.cam.x$mantel.res[,3], col=7, lty=3)
for(i in 1:length(mc.cam.x$mantel.res[,1])){
  if(is.na(mc.cam.x$mantel.res[i,5]==TRUE)){print("NA")
  }else{if(mc.cam.x$mantel.res[i,5]<.05){points(mc.cam.x$mantel.res[i,1],mc.cam.x$mantel.res[i,3],pch=22,col=7,bg=7)}}} 
    