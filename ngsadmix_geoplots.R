#!/bin/sh
#!/usr/bin/env Rscript
## module load system/R-3.6.1 ; R
# module load system/R-3.5.2 ; module load system/R-3.5.2 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R
module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R
###module load system/R-3.5.2 ; module load system/R-3.5.2 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R

##to be deleted when into script
awk '{print $6
for loc in xchrom ; do #autosomes
  prefix=rad.ptri
  order_file=$data/pango_sample_list.txt
  sample_file=pango.rad.all_$loc
  ngsadmix_dir=$res/ngsadmix.$loc
  module load system/R-3.4.3
  Rscript --vanilla $bin/ngsadmix_plots.R $popmap2 $prefix.${sample_file} $ngsadmix_dir $order_file
  module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3
  Rscript --vanilla $bin/ngsadmix_geoplots.R $popmap2 $prefix.${sample_file} $ngsadmix_dir $order_file
done

########popmap=popmap2

ngsadmix_dir="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom"
sample_file="pango.rad.all_xchrom"
order_file="/work/idumville/pangolins/RAD/data/pango_sample_list.txt"
prefix_sample_file="rad.ptri.pango.rad.all_xchrom"
prefix="rad.ptri"
popmap="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom/sampledataxchrom.txt"
loc="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom"


#for autosomes
ngsadmix_dir="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes"
sample_file="pango.rad.all_autosomes"
order_file="/work/idumville/pangolins/RAD/data/pango_sample_list.txt"
prefix_sample_file="rad.ptri.pango.rad.all_autosomes"
prefix="rad.ptri"
popmap="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/sampledataautosomes.txt"
loc="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes"

#for lowread autosomes
ngsadmix_dir="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads"
sample_file="pango.rad.all_autosomes"
order_file="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads/highreadsamples.txt"
prefix_sample_file="lowreads.pango.rad.autosomes"
prefix="lowreads"
popmap="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads/sampledataautosomes.txt"
loc="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads"


#for xchrom lowreads
ngsadmix_dir="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom.lowreads"
sample_file="pango.rad.all_xchrom"
order_file="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads/highreadsamples.txt"
prefix_sample_file="lowreads.pango.rad.xchrom"
prefix="lowreads"
popmap="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom.lowreads/sampledataxchrom.txt"
loc="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom.lowreads"



##end of to be deleted when into script

args = commandArgs(trailingOnly=TRUE)
prefix_sample_file<-args[2]
loc<-args[3]
order_file<-args[4]
# prefix_sample_file<-"Xrad.mms.rad.all" ; args[1]<-"../../data/RAD_PG2020.data.info" ; loc="/work/pgaubert/monk_seals/RAD/results/ngsadmix/" ;order_file<-"../../data/RAD_PG2020_sample_list.txt"
setwd(loc)


.libPaths( "/tools/R/R-3.4.3/lib64/R/library")


# call libraries 
library(raster)
library(rgdal)
library(ggplot2)
library(reshape2)
library(plyr)
library(scatterpie)
library(dplyr)
library(rgdal)
library(raster) 
library(rworldmap) 
library(ggrepel)
library(RColorBrewer)
library(raster)
library(ggnewscale)
library(png)
library(cowplot)




# call likevalues
likevalues<-read.table(paste("./",prefix_sample_file,"_likevalues.txt",sep=""))

psmall <- readRDS("../../data/geodata/smallmap.rdata")
pbig <- readRDS("../../data/geodata/bigmap.rdata")
psmall_np <- readRDS("../../data/geodata/smallmap_noplane.rdata")
pbig_np <- readRDS("../../data/geodata/bigmap_noplane.rdata")

##########################################
####### Plotting rural  ########## 
##########################################
##extract out of here for final plot

pdf(paste("./",prefix_sample_file,"_geo_rural.pdf",sep=""), 6, 6) #10,6
#par(mai = c(2, 2, 2, 5))
for(j in 2:12){
  samples<-read.table(popmap,header=F)
  #for X only 
   #samples<-samples[-which(samples$V1 == "Y199_1"),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
   urbsamples<-which(samples$V10 %in% c("Urban_Market", "Seizure")) # %in% c("Yaounde", "Abgab", "Tongue", "Benin", "DCM", "Cote_dIvoire", "DDM", "Ghana", "SeizureP", "SeizureB", "Central_African_Republic"))
   westsamples<-which(samples$V2 %in% c("Kumasi", "Mankessim"))
   samples<-samples[-c(urbsamples),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
   samples<-samples[-c(westsamples),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
  popx<-samples[,c(8:7,2)] ; popfor<-aggregate(popx[,1:2], list(popx[,3]), mean) #long:lat, popn
bestlike<-which(likevalues[likevalues[,2]==j,1]==max(likevalues[likevalues[,2]==j,1]))
admix<-as.matrix(read.table(paste("./",prefix_sample_file,"_k",j,"/",prefix_sample_file,"_k",j,"_s",bestlike,".qopt",sep=""))) 
rownames(admix)<-seq(length=nrow(admix)) ;
admix<-admix[-c(urbsamples),]
admix<-admix[-c(westsamples),]
admix<-aggregate(admix[,1:ncol(admix)], list(samples[,2]), sum); 
admix[,ncol(admix)+1]<-NA#admix[,ncol(admix)]
for(i in 1:length(admix[,1])){ admix[i,ncol(admix)]<-sum(admix[i,2:(ncol(admix)-1)]) }
popfor[,4:(ncol(popfor)+(ncol(admix)-1))]<-admix[,2:ncol(admix)]
colnames(popfor)[4:ncol(popfor)]<-c(LETTERS[1:(ncol(admix)-2)],"radius")

#plot 

#set smaller limits
lon <- c(8.000, 14.2194) #-10.5318, 29.53505) #
lat <- c(-0.271891, 7.25) #-3.7346, 10.7618) #

p2 <- psmall_np + new_scale_color() + new_scale_fill() + scale_fill_brewer(palette = "Set3") + scale_color_brewer() + geom_scatterpie(data=popfor, aes(x=V8, y=V7, r=.02*radius), cols=LETTERS[1:(ncol(admix)-2)], alpha=1, color="black")+ geom_scatterpie(data=popfor, aes(x=V8, y=V7, r=.0199*radius), cols=LETTERS[1:(ncol(admix)-2)], color=NA) 

#getting lineage informaiton onto dataframe
popfor <- cbind(popfor, lineage=samples[match(popfor[,1], samples[,2]),9])
p2 <- p2 + geom_text_repel(data=popfor, aes(x=V8, y=V7, label=Group.1), min.segment.length = 0.1, size=2, segment.size=0.1, fontface="bold", bg.color = "white", bg.r = .01)
p2 <- p2  + coord_fixed() +
  theme(legend.position = "null") + ggtitle(paste0("Rural; k=", j)) 

print(p2) #1*radius
}

dev.off()

###################
### Final WCA plot ###
###################3
 #10,6
#par(mai = c(2, 2, 2, 5))
j=6
  samples<-read.table(popmap,header=F)
   samples$V2 <- as.character(samples$V2)
  samples$V2[which(samples$V2 == "SeizureR")] <- "SeizureP"
  samples$V2[which(samples$V2 == "Gnkoltang")] <- "Nkoltang"
  samples$V2[which(samples$V2 == "GMedouneu")] <- "Medouneu"
 westsamples<-which(samples$V2 %in% c("Kumasi", "Mankessim"))
   samples<-samples[-c(which(samples$V2 %in% c("Kumasi", "Mankessim"))),]
   samples$V2 <- as.factor(samples$V2)
   #samples<-droplevels(samples) ; rownames(samples) <- NULL
  popx<-samples[,c(8:7,2)] ; popfor<-aggregate(popx[,1:2], list(popx[,3]), mean) #long:lat, popn
bestlike<-which(likevalues[likevalues[,2]==j,1]==max(likevalues[likevalues[,2]==j,1]))
admix<-as.matrix(read.table(paste("./",prefix_sample_file,"_k",j,"/",prefix_sample_file,"_k",j,"_s",bestlike,".qopt",sep=""))) 
rownames(admix)<-seq(length=nrow(admix)) ;

admix<-admix[-c(westsamples),]
admix<-aggregate(admix[,1:ncol(admix)], list(samples[,2]), sum); 
admix[,ncol(admix)+1]<-NA#admix[,ncol(admix)]
for(i in 1:length(admix[,1])){ admix[i,ncol(admix)]<-sum(admix[i,2:(ncol(admix)-1)]) }
popfor[,4:(ncol(popfor)+(ncol(admix)-1))]<-admix[,2:ncol(admix)]
colnames(popfor)[4:ncol(popfor)]<-c(LETTERS[1:(ncol(admix)-2)],"radius")
#getting lineage informaiton onto dataframe
popfor <- cbind(popfor, lineage=samples[match(popfor[,1], samples[,2]),9])

popfor <- popfor[-which(popfor$Group.1 %in% c("Yaounde_Nkolndongo_Market", "Douala_Dakat_Market", "SeizureR", "SeizureP", "Douala_Central_Market", "SeizureB", "Banyo")),]

#plot 
lon <- c(8.000, 14.2194) 
lat <- c(-0.271891, 7.25) 


png(paste("./Figure3A",prefix_sample_file,"_geo_all_noscale_final.png",sep=""),units="in", 2,2.4, res=900) 
p2 <- psmall_np + new_scale_color() + new_scale_fill() + scale_fill_brewer(palette = "Set3") + scale_color_brewer() +  geom_scatterpie(data=popfor, aes(x=V8, y=V7, r=.20), cols=LETTERS[1:(ncol(admix)-2)], alpha=1, color="black")+ geom_scatterpie(data=popfor, aes(x=V8, y=V7, r=.199), cols=LETTERS[1:(ncol(admix)-2)], color=NA)
p2 <- p2 + geom_text_repel(data=popfor, aes(x=V8, y=V7, label=Group.1), min.segment.length = 0.1, size=1.2, segment.size=0.1, fontface="bold", bg.color = "white", bg.r = .05, box.padding = 0.05)
p2 <- p2  + coord_fixed() +
  theme(legend.position = "null", axis.text=element_text(size=5), plot.margin = unit(c(-0.5, 0.1, -1, -0.4), 'cm')) + 
  annotate(geom="text", x=8.25, y=7, label="A", colour='black',  fontface = "bold")   
print(p2)
dev.off()

#, plot.margin = unit(c(0, 0.5, -1, -0.4), 'cm')

low_list="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/lowurbanreads.txt"

png(paste("./Figure3B",prefix_sample_file,"_ngsadmix_rural_barplot_final.png",sep=""), units="in", 0.6, 2.4, res=900) #60,30 for pdf in inch; 6000,3000 pixels png
layout(matrix(c(1,1), nrow = 1, ncol = 1, byrow = TRUE))
i=6
K=i
  samples<-read.table(popmap,header=F)
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
      poplist <- c("SeizureP", "SeizureB", "Gnkoltang", "Oyem", "Makokou", "Medouneu", "Bisobinam", "Taguete","Misergue", "Emangos", "Anguma", "Ebenguan", "Bongoro", "Campo", "Maan", "Djoum", "Sangmelima", "Lolodorf", "Bipindi", "Ekombitie", "Eseka", "Akonolinga", "Yaounde_Nkolndongo_Market", "Esse", "Abong_Mbang",  "Douala_Central_Market", "Douala_Dakat_Market", "Yabassi", "Bayomen", "Manengole", "Foumbot", "Bayib-Assibong", "Nditam", "Manyemen_Nguti",  "Banyo","Bioko_Batete", "Mankessim",  "Kumasi") #
  samples[,2]<-factor(samples[,2],levels=poplist) 
  admix<-admix[,order(samples[,2])]; 
  samples<-samples[order(samples[,2]),]
  niveaux<-levels(as.factor(samples[,2]))
  pop<-samples[order(match(samples[,2],niveaux)),]
  
  #HERE REMOVE INDIVIDUALS FROM POP, SAMPLES AND ADMIX FOR RURAL HIGH READS 
  lines <- readLines(low_list)
  lines <- c(lines, "By3_1")
  samplenum <- which(samples$V1 %in% lines)
  samples <- samples[-c(samplenum),]
  pop <- pop[-c(samplenum),]
  admix <- admix[,-c(samplenum)]

  tempo<-tapply(1:nrow(pop),pop[,2],mean);
  tempo<-tempo[order(tempo)]
  tempomax<-tapply(1:nrow(pop),pop[,2],max);
  tempomax<-tempomax[order(tempomax)]
  
    par(mar=c(0.1,0.1,0,0.1))
      barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samples$V1)
      abline(h=c(0,tempomax),lty=1,lwd=0.6,col="grey50")
    text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,col="black", cex=0.2) #cex=1.5,
      
      mtext("B",side=3,line=-0.5,adj=0, padj=1, font = 2)  
dev.off()

#trying to get both on one plot, still a WIP
#library(grid) ; library(ggplotify)

vp.BottomRight <- viewport(#height=5, width=3.5,  default.units = "inches",
                           just=1, 
                           y=0.5, x=0.5)

par(mfrow=c(1,2))
png("./test.png")
par(omi=c(0.1,1,0.1,2.5))
barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samples$V1, cex.names=0.1)  #names.arg=samples$coverage #, cex.names=1 
    text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,col="black", cex=0.2) #cex=1.5,
      abline(h=c(0,tempomax),lty=1,lwd=1,col="grey22")  
      p3 <- recordPlot()

# plot the ggplot using the print command
print(p2, vp=vp.BottomRight)

png("./test.png")
par(mfrow=c(1,2))
print(p3, p2, ncol=2)
dev.off()

#need ggplotify
as.ggplot()

##########################################
####### Plotting all  ########## 
##########################################

pdf(paste("./",prefix_sample_file,"_geo_all_noscale.pdf",sep=""), 6, 6)

#par(mai = c(2, 2, 2, 5))
for(j in 2:12){ #colour palette only to 12 
  samples<-read.table(popmap,header=F)
  #for X only 
   #samples<-samples[-which(samples$V1 == "Y199_1"),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
   
   westsamples<-which(samples$V2 %in% c("Kumasi", "Mankessim"))
   
   samples<-samples[-c(westsamples),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
  popx<-samples[,c(8:7,2)] ; popfor<-aggregate(popx[,1:2], list(popx[,3]), mean) #long:lat, popn
bestlike<-which(likevalues[likevalues[,2]==j,1]==max(likevalues[likevalues[,2]==j,1]))
admix<-as.matrix(read.table(paste("./",prefix_sample_file,"_k",j,"/",prefix_sample_file,"_k",j,"_s",bestlike,".qopt",sep=""))) 
rownames(admix)<-seq(length=nrow(admix)) ;

admix<-admix[-c(westsamples),]
admix<-aggregate(admix[,1:ncol(admix)], list(samples[,2]), sum); 
admix[,ncol(admix)+1]<-NA#admix[,ncol(admix)]
for(i in 1:length(admix[,1])){ admix[i,ncol(admix)]<-sum(admix[i,2:(ncol(admix)-1)]) }
popfor[,4:(ncol(popfor)+(ncol(admix)-1))]<-admix[,2:ncol(admix)]
colnames(popfor)[4:ncol(popfor)]<-c(LETTERS[1:(ncol(admix)-2)],"radius")

#plot 


#set smaller limits
lon <- c(8.000, 14.2194) #-10.5318, 29.53505) #
lat <- c(-0.271891, 7.25) #-3.7346, 10.7618) #

p2 <- psmall + new_scale_color() + new_scale_fill() + scale_fill_brewer(palette = "Set3") + scale_color_brewer() +  geom_scatterpie(data=popfor, aes(x=V8, y=V7, r=.20), cols=LETTERS[1:(ncol(admix)-2)], alpha=1, color="black")+ geom_scatterpie(data=popfor, aes(x=V8, y=V7, r=.199), cols=LETTERS[1:(ncol(admix)-2)], color=NA)

#getting lineage informaiton onto dataframe
popfor <- cbind(popfor, lineage=samples[match(popfor[,1], samples[,2]),9])
p2 <- p2 + geom_text_repel(data=popfor, aes(x=V8, y=V7, label=Group.1), min.segment.length = 0.1, size=2, segment.size=0.1, fontface="bold", bg.color = "white", bg.r = .01)
p2 <- p2  + coord_fixed() +
  theme(legend.position = "null") + ggtitle(paste0("All; k=", j)) 
 
print(p2) #1*radius
}

  
dev.off()




##########################################
####### Figure 2 AL autosomes  ########## 
##########################################

ngsadmix_dir="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes"
sample_file="pango.rad.all_autosomes"
order_file="/work/idumville/pangolins/RAD/data/pango_sample_list.txt"
prefix_sample_file="rad.ptri.pango.rad.all_autosomes"
prefix="rad.ptri"
popmap="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/sampledataautosomes.txt"
loc="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes"
setwd(loc)

lon <- c(-10.5318, 29.53505) #
lat <- c(-3.7346, 10.7618) #


#Loading PCA
PCAall<-readPNG("../pcangsd.autosomes/rad.ptri.pango.rad.all_autosomes.allpcangsd_first_plots2.png")

likevalues<-read.table(paste("./",prefix_sample_file,"_likevalues.txt",sep=""))


  samples<-read.table(popmap,header=F)
  samples$V2 <- as.character(samples$V2)
  samples$V2[which(samples$V2 == "SeizureR")] <- "SeizureP"
  samples$V2[which(samples$V2 == "Gnkoltang")] <- "Nkoltang"
  samples$V2[which(samples$V2 == "GMedouneu")] <- "Medouneu"
  samples$V2 <- as.factor(samples$V2)
  
popfor <- popfor[-which(popfor$Group.1 %in% c("Yaounde_Nkolndongo_Market", "Douala_Dakat_Market", "SeizureR", "SeizureP", "Douala_Central_Market", "SeizureB", "Banyo")),]

png("./Figure2A.png",units="in", 6.5, 2.65, res=900)
j=6
  popx<-samples[,c(8:7,2)] ; popfor<-aggregate(popx[,1:2], list(popx[,3]), mean) #long:lat, popn
bestlike<-which(likevalues[likevalues[,2]==j,1]==max(likevalues[likevalues[,2]==j,1]))
admix<-as.matrix(read.table(paste("./",prefix_sample_file,"_k",j,"/",prefix_sample_file,"_k",j,"_s",bestlike,".qopt",sep=""))) 
rownames(admix)<-seq(length=nrow(admix)) ;
#exporting the likelihoods per sample -> 
#rownames(admix) <- samples$V1
#write.table(admix, file="./ALlikelihood.txt", quote=FALSE, row.names=FALSE
admix<-aggregate(admix[,1:ncol(admix)], list(samples[,2]), sum); 
admix[,ncol(admix)+1]<-NA #admix[,ncol(admix)]
for(i in 1:length(admix[,1])){ admix[i,ncol(admix)]<-sum(admix[i,2:(ncol(admix)-1)]) }
popfor[,4:(ncol(popfor)+(ncol(admix)-1))]<-admix[,2:ncol(admix)]
colnames(popfor)[4:ncol(popfor)]<-c(LETTERS[1:(ncol(admix)-2)],"radius")
#getting lineage informaiton onto dataframe 
popfor <- cbind(popfor, lineage=samples[match(popfor[,1], samples[,2]),9])

popfor <- popfor[-which(popfor$Group.1 %in% c("Banyo", "Douala_Central_Market", "Douala_Dakat_Market", "Yaounde_Nkolndongo_Market")),]

p2 <- pbig_np  + new_scale_color() + new_scale_fill() + scale_fill_brewer(palette = "Set3") + scale_color_brewer() +  geom_scatterpie(data=popfor, aes(x=V8, y=V7, r=.20), cols=LETTERS[1:(ncol(admix)-2)], alpha=1, color="black")+ geom_scatterpie(data=popfor, aes(x=V8, y=V7, r=.199), cols=LETTERS[1:(ncol(admix)-2)], color=NA)

p2 <- p2 + new_scale_color()  + scale_color_brewer(palette="Set2") #+ geom_text_repel(data=popfor, aes(x=V8, y=V7, label=Group.1, color = lineage), box.padding=0.05, min.segment.length = 0.1, size=1.2, segment.size=0.1, fontface="bold", bg.color = "black", bg.r = .1, max.overlaps=Inf)
p2 <- p2  + coord_fixed() +
  theme(legend.position = "null",plot.margin = unit(c(0,0.1,-1,-0.3), 'cm'))  +  annotate(geom="text", x=-10, y=10, label="A", colour='black',  fontface = "bold") 
plot(p2)
dev.off()


##########################################
####### Plotting xchrom  ########## 
##########################################
ngsadmix_dir="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom"
sample_file="pango.rad.all_xchrom"
order_file="/work/idumville/pangolins/RAD/data/pango_sample_list.txt"
prefix_sample_file="rad.ptri.pango.rad.all_xchrom"
prefix="rad.ptri"
popmap="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom/sampledataxchrom.txt"
loc="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom"

lon <- c(-10.5318, 29.53505) #
lat <- c(-3.7346, 10.7618) #


likevalues<-read.table(paste("./",prefix_sample_file,"_likevalues.txt",sep=""))

pdf(paste("./",prefix_sample_file,"_geo_all_noscale.pdf",sep=""), 10, 6)  #10,6

for(j in 2:10){
  samples<-read.table(popmap,header=F)
  samples$V2[which(samples$V2 == "SeizureR")] <- "SeizureP"
   samples<-samples[-which(samples$V1 == "Y199_1"),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
  popx<-samples[,c(8:7,2)] ; popfor<-aggregate(popx[,1:2], list(popx[,3]), mean) #long:lat, popn
bestlike<-which(likevalues[likevalues[,2]==j,1]==max(likevalues[likevalues[,2]==j,1]))
admix<-as.matrix(read.table(paste("./",prefix_sample_file,"_k",j,"/",prefix_sample_file,"_k",j,"_s",bestlike,".qopt",sep=""))) 
rownames(admix)<-seq(length=nrow(admix)) ;

admix<-aggregate(admix[,1:ncol(admix)], list(samples[,2]), sum); 
admix[,ncol(admix)+1]<-NA#admix[,ncol(admix)]
for(i in 1:length(admix[,1])){ admix[i,ncol(admix)]<-sum(admix[i,2:(ncol(admix)-1)]) }
popfor[,4:(ncol(popfor)+(ncol(admix)-1))]<-admix[,2:ncol(admix)]
colnames(popfor)[4:ncol(popfor)]<-c(LETTERS[1:(ncol(admix)-2)],"radius")
#popfor$radius <- (popfor$radius)
#plot 

p2 <- pbig + new_scale_color() + new_scale_fill() + scale_fill_brewer(palette = "Set3") + scale_color_brewer() +  geom_scatterpie(data=popfor, aes(x=V8, y=V7, r=.20), cols=LETTERS[1:(ncol(admix)-2)], alpha=1, color="black")+ geom_scatterpie(data=popfor, aes(x=V8, y=V7, r=.199), cols=LETTERS[1:(ncol(admix)-2)], color=NA) 
#getting lineage informaiton onto dataframe
popfor <- cbind(popfor, lineage=samples[match(popfor[,1], samples[,2]),9])
p2 <- p2 + new_scale_color()  + scale_color_brewer(palette="Set2") + geom_text_repel(data=popfor, aes(x=V8, y=V7, label=Group.1, color = lineage), min.segment.length = 0.1, size=1, segment.size=0.1, fontface="bold", bg.color = "black", bg.r = .1, max.overlaps=Inf)
p2 <- p2  + coord_fixed() +
  theme(legend.position = "null") + ggtitle(paste0("Xchrom; k=", j)) 
print(p2) #1*radius
}

  
dev.off()





exit






##########################################
####### scrap ########## 
##########################################

#ancient <- shapefile("/work/pgaubert/monk_seals/RAD/data/gis/Ancient_distribution_M._monachus1-polygon.shp")
#ancient_crs <- spTransform(ancient,CRS("+proj=longlat +datum=WGS84")) #change GCS
#ancient_df <- fortify(ancient_crs) #fortify shpfile
#save(ancient_df,file="ancient.rda")
#load("./ancient.rda")
#ancient_clip <- raster::intersect(ancient, clipper_med)
#ancient_clip_f <- fortify(ancient_clip)
  #  samples<-samples[-which(samples$V1 == "Y199_1"),]
#  samples<-droplevels(samples)
#  rownames(samples) <- NULL
  #samples<-read.table(args[1],header=T); sort(samples[,3])
  #pg2020<-read.table(order_file,header=F)  
  #samples[,2]<-factor(samples[,2],levels=poplist);
  
  #admix<-admix[,order(pg2020)]; #[,1]
  
  #samples<-samples[order(samples[,2]),]
  
  #admix<-admix[,order(samples[,2])] #[,2]
  
  #admix<-admix[,order(rownames(samples))]
  
  #admix<-admix[,samples$V5>2] #coverage
  #samples<-samples[samples$V5>2,] #coverage
#adding labels of places TO EDIT
seiz <- grobTree(textGrob("Seizures", x=0.165,  y=0.22, hjust=0,gp=gpar(col="black", fontsize=11)))
mad <- grobTree(textGrob("Madeira", x=0.17,  y=0.55, hjust=0,gp=gpar(col=rainbow(4)[4], fontsize=13)))
wmed <- grobTree(textGrob("W-Med", x=0.47,  y=0.7, hjust=0,gp=gpar(col="darkgreen", fontsize=13)))
cmed <- grobTree(textGrob("C-Med", x=0.66,  y=0.57, hjust=0,gp=gpar(col=rainbow(4)[1], fontsize=13)))
emed <- grobTree(textGrob("E-Med", x=0.78,  y=0.52, hjust=0,gp=gpar(col=rainbow(4)[1], fontsize=13)))
bsea <- grobTree(textGrob("Black-Sea", x=0.85,  y=0.79, hjust=0,gp=gpar(col="darkgreen", fontsize=13)))

 #+ 
  #      annotation_custom(seiz)+ 
   #     annotation_custom(mad)+
    #    annotation_custom(wmed)+
     #   annotation_custom(cmed)+
      #  annotation_custom(emed)+
       # annotation_custom(bsea))
       
       
       
###### Represent clustering results geographically
# Mapping settings
world <- getMap(resolution = "high")
lon <- c(-12, 30)
lat <- c(-4, 11) #chosen for seizures: -0.599851	3.032364 (middle of ocean) used to be yaounde urban market coords
clipper_med <- as(extent(lon, lat), "SpatialPolygons")
proj4string(clipper_med) <- CRS(proj4string(world))
world_clip <- raster::intersect(world, clipper_med)
world_clip_f <- fortify(world_clip)


pdf(paste("./",prefix_sample_file,"_ngsadmix_geoplot.pdf",sep=""), 10, 6) #_names
for(K in 2:8){
K=6 ; 
i=K
  samples<-read.table(popmap,header=F); sort(samples[,1])
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
 admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
 #for k=6, saving admix to plo; rotating and adding names
  poplist<-c("Guinea", "Cote_dIvoire", "Ghana", "Togo", "Benin", "Nigeria", "Banyo", "Foumbot", "Nditam", "Manengole", "SO", "Yabassi", "DCM", "DDM", "Bioko", "Eseka", "Bayomen", "Yaounde", "Tongue", "Ekombitie", "Akonolinga", "Abgab", "Abongmbang", "Sangmelima",  "PNCM", "GES", "NorthG", "GFV", "Central_African_Republic", "DR_Congo", "SeizureB", "SeizureP")


MMS<-samples[,c(7,8)] ; colnames(MMS)[1]<-"lati" ; colnames(MMS)[2]<-"long"
MMS[,1]<-as.numeric(MMS[,1])
MMS[,2]<-as.numeric(MMS[,2])
MMS[,3]<-samples[,2] ; colnames(MMS)[3]<-"pop"
MMS[,4:(K+3)]<-t(admix);colnames(MMS)[4:(K+3)]<-LETTERS[1:K]

MMS$lineage<-samples[,9]
MMS<-na.omit(MMS,na.strings = c(NA, "NA"))
library(dplyr)
MMS <- add_count(MMS, pop)
MMS_map <- ggplot() + 
	geom_rect(mapping=aes(xmin=lon[1], xmax=lon[2], ymin=lat[1], ymax=lat[2]), color="lightblue",fill="lightblue", alpha=0.5)+

  geom_polygon(data = world_clip_f, aes(x = long, y = lat, group = group),
		fill = "lightgrey", colour = "grey") + 
  xlim(lon[1], lon[2]) + ylim(lat[1], lat[2])+
  coord_quickmap()+
	theme_bw() +
  geom_scatterpie(aes(x=long, y=lati), data=MMS, cols=c(LETTERS[1:K]),color=NA,alpha=.7)+ 
  
  geom_text_repel(data=MMS[match(unique(MMS[,2]), MMS[,2]),], aes(x=long, y=lati, label=pop, colour=lineage), point.size = NA, min.segment.length = 0, size=2) + 
	xlab("") +
	ylab("") +
  theme(legend.position = "right",plot.margin = unit(c(0, 0, -1, -0.8), 'lines'),
        panel.border = element_rect(colour = 'black', fill = 'transparent'))+
	xlab("") +
 scale_color_brewer(palette = "Dark2") +
	ylab("")  
library(grid)

print(MMS_map)

}
dev.off()

#plotting pie as population

pdf(paste("./",prefix_sample_file,"_ngsadmix_geoplot_pop.pdf",sep=""), 15, 15)
  
for(j in 2:10){
  samples<-read.table(popmap,header=F); sort(samples[,1])
  #for X only 
   samples<-samples[-which(samples$V1 == "Y199_1"),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
  popx<-samples[,c(8:7,2)] ; popfor<-aggregate(popx[,1:2], list(popx[,3]), mean) #long:lat, popn
bestlike<-which(likevalues[likevalues[,2]==j,1]==max(likevalues[likevalues[,2]==j,1]))
admix<-as.matrix(read.table(paste("./",prefix_sample_file,"_k",j,"/",prefix_sample_file,"_k",j,"_s",bestlike,".qopt",sep=""))) 
rownames(admix)<-seq(length=nrow(admix)) ;
admix<-aggregate(admix[,1:ncol(admix)], list(samples[,2]), sum); 
admix[,ncol(admix)+1]<-NA#admix[,ncol(admix)]
for(i in 1:length(admix[,1])){ admix[i,ncol(admix)]<-sum(admix[i,2:(ncol(admix)-1)]) }
popfor[,4:(ncol(popfor)+(ncol(admix)-1))]<-admix[,2:ncol(admix)]
colnames(popfor)[4:ncol(popfor)]<-c(LETTERS[1:(ncol(admix)-2)],"radius") ;
#plot 
p2 <- ggplot() + 
	geom_rect(mapping=aes(xmin=lon[1], xmax=lon[2], ymin=lat[1], ymax=lat[2]), color="lightblue",fill="lightblue", alpha=0.5)+
  geom_polygon(data = world_clip_f, aes(x = long, y = lat, group = group),
		fill = "lightgrey", colour = "grey") + 
  xlim(lon[1], lon[2]) + ylim(lat[1], lat[2])+
  coord_quickmap()+
	theme_bw()
p2 <- p2 + geom_scatterpie(aes(x=V8, y=V7, r=.003*radius), data=popfor, cols=LETTERS[1:(ncol(admix)-2)], color=NA, alpha=.7)
p2 <- p2 + geom_text_repel(data = popfor, aes(x = V8, y = V7, label = Group.1), size = 3, point.padding=0.4,colour="yellow", fontface="bold")
p2 <- p2
#p2 + geom_scatterpie_legend(c(.015,.030,.060), x=49, y=-16.6,n=3)
print(p2) #1*radius

#p <- grid.arrange(p1, p2, nrow = 1)
#print(p)
}
dev.off()


#another plot without yaounde
pdf(paste("./",prefix_sample_file,"_ngsadmix_geoplot_pop_woYaounde.pdf",sep=""), 10, 6)
  
for(j in 2:10){
  samples<-read.table(popmap,header=F); sort(samples[,1])
  #for X only 
   #samples<-samples[-which(samples$V1 == "Y199_1"),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
   yaosamples<-which(samples$V2 == "Yaounde")
   samples<-samples[-which(samples$V2 == "Yaounde"),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
  popx<-samples[,c(8:7,2)] ; popfor<-aggregate(popx[,1:2], list(popx[,3]), mean) #long:lat, popn
bestlike<-which(likevalues[likevalues[,2]==j,1]==max(likevalues[likevalues[,2]==j,1]))
admix<-as.matrix(read.table(paste("./",prefix_sample_file,"_k",j,"/",prefix_sample_file,"_k",j,"_s",bestlike,".qopt",sep=""))) 
rownames(admix)<-seq(length=nrow(admix)) ;
admix<-admix[-c(yaosamples),]
admix<-aggregate(admix[,1:ncol(admix)], list(samples[,2]), sum); 
admix[,ncol(admix)+1]<-NA#admix[,ncol(admix)]
for(i in 1:length(admix[,1])){ admix[i,ncol(admix)]<-sum(admix[i,2:(ncol(admix)-1)]) }
popfor[,4:(ncol(popfor)+(ncol(admix)-1))]<-admix[,2:ncol(admix)]
colnames(popfor)[4:ncol(popfor)]<-c(LETTERS[1:(ncol(admix)-2)],"radius") ;
popfor$radius <- (popfor$radius + 20)*2
#plot 
p2 <- ggplot() + 
	geom_rect(mapping=aes(xmin=lon[1], xmax=lon[2], ymin=lat[1], ymax=lat[2]), color="lightblue",fill="lightblue", alpha=0.5)+
  geom_polygon(data = world_clip_f, aes(x = long, y = lat, group = group),
		fill = "lightgrey", colour = "grey") + 
  xlim(lon[1], lon[2]) + ylim(lat[1], lat[2])+
  coord_quickmap()+
	theme_bw()
p2 <- p2 + geom_scatterpie(aes(x=V8, y=V7, r=.003*radius), data=popfor, cols=LETTERS[1:(ncol(admix)-2)], color=NA, alpha=.7)
#getting lineage informaiton onto dataframe
popfor <- cbind(popfor, lineage=samples[match(popfor[,1], samples[,2]),9])
p2 <- p2 + geom_text_repel(data=popfor, aes(x=V8, y=V7, label=Group.1, colour=lineage), min.segment.length = 0.5, size=2, segment.size=0.1) # geom_text_repel(data = popfor, aes(x = V8, y = V7, label = Group.1), size=2, colour="yellow") #
p2 <- p2 + 
  theme(legend.position = "bottom",plot.margin = unit(c(0, 0, -1, -0.8), 'lines'),
        panel.border = element_rect(colour = 'black', fill = 'transparent')) + ggtitle(paste0("Without Yaounde; k=", j))
#p2 + geom_scatterpie_legend(c(.015,.030,.060), x=49, y=-16.6,n=3)
print(p2) #1*radius

#p <- grid.arrange(p1, p2, nrow = 1)
#print(p)
}
dev.off()



##########################################
####### Plotting urban without radius scaling ########## 
##########################################
pdf(paste("./",prefix_sample_file,"_geo_urban_noscale.pdf",sep=""), 10, 6)  #10,6
#par(mai = c(2, 2, 2, 5))
for(j in 2:10){
  samples<-read.table(popmap,header=F)
  #for X only 
   #samples<-samples[-which(samples$V1 == "Y199_1"),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
   urbsamples<-which(samples$V10 %in% c("Urban_Market", "Seizure")) # %in% c("Yaounde", "Abgab", "Tongue", "Benin", "DCM", "Cote_dIvoire", "DDM", "Ghana", "SeizureP", "SeizureB", "Central_African_Republic"))
   samples<-samples[c(urbsamples),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
  popx<-samples[,c(8:7,2)] ; popfor<-aggregate(popx[,1:2], list(popx[,3]), mean) #long:lat, popn
bestlike<-which(likevalues[likevalues[,2]==j,1]==max(likevalues[likevalues[,2]==j,1]))
admix<-as.matrix(read.table(paste("./",prefix_sample_file,"_k",j,"/",prefix_sample_file,"_k",j,"_s",bestlike,".qopt",sep=""))) 
rownames(admix)<-seq(length=nrow(admix)) ;
admix<-admix[c(urbsamples),]
admix<-aggregate(admix[,1:ncol(admix)], list(samples[,2]), sum); 
admix[,ncol(admix)+1]<-NA#admix[,ncol(admix)]
for(i in 1:length(admix[,1])){ admix[i,ncol(admix)]<-sum(admix[i,2:(ncol(admix)-1)]) }
popfor[,4:(ncol(popfor)+(ncol(admix)-1))]<-admix[,2:ncol(admix)]
colnames(popfor)[4:ncol(popfor)]<-c(LETTERS[1:(ncol(admix)-2)],"radius")
#popfor$radius <- (popfor$radius)
#plot 
p2 <- ggplot() + #data = rast.table, aes(x = x, y = y)
  #geom_tile(fill = rast.table$rgb) +
  #scale_alpha_discrete(range=c(1,0)) +
  geom_path(data=count.subset, aes(x = long, y = lat, group = group), color = 'black', size=0.5) +
  geom_path(data=river.subset, aes(x = long, y = lat, group = group), color = 'blue', alpha=0.7, size=0.3) +
  geom_path(data=road.subset, aes(x = long, y = lat, group = group), color = 'azure3', alpha=0.7, size=0.3) +
  geom_path(data=rail.subset, aes(x = long, y = lat, group = group), color = 'azure2', alpha=0.7, size=0.3) +
  geom_polygon(data=urban.subset, aes(x = long, y = lat, group = group), fill = "beige", colour = "beige") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('') + ylab('')
p2 <- p2 + geom_scatterpie(aes(x=V8, y=V7, r=0.25), data=popfor, cols=LETTERS[1:(ncol(admix)-2)], alpha=1, color=NA) #, r=.015*radius no radius on this one
#getting lineage informaiton onto dataframe
popfor <- cbind(popfor, lineage=samples[match(popfor[,1], samples[,2]),9])
p2 <- p2 + geom_text_repel(data=popfor, aes(x=V8, y=V7, label=Group.1, colour=lineage), min.segment.length = 0.1, size=1, segment.size=0.1, bg.color = "black", bg.r = .05) # geom_text_repel(data = popfor, aes(x = V8, y = V7, label = Group.1), size=2, colour="yellow") #, fontface="bold"
p2 <- p2 + theme_bw() +
  theme(legend.position = "bottom",plot.margin = unit(c(0, 0.5, 0, 0), 'lines'),
        panel.border = element_rect(colour = 'black', fill = 'transparent')) + ggtitle(paste0("Without Rural; k=", j))
#p2 + geom_scatterpie_legend(c(.015,.030,.060), x=49, y=-16.6,n=3)
print(p2) #1*radius
}

  
dev.off()


##########################################
####### Plotting rural without radius scaling ########## 
##########################################
pdf(paste("./",prefix_sample_file,"_geo_rural_noscale.pdf",sep=""),6, 6) #10,6
#par(mai = c(2, 2, 2, 5))
for(j in 2:10){
  samples<-read.table(popmap,header=F)
  #for X only 
   #samples<-samples[-which(samples$V1 == "Y199_1"),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
   urbsamples<-which(samples$V10 %in% c("Urban_Market", "Seizure")) # %in% c("Yaounde", "Abgab", "Tongue", "Benin", "DCM", "Cote_dIvoire", "DDM", "Ghana", "SeizureP", "SeizureB", "Central_African_Republic"))
   samples<-samples[-c(urbsamples),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
  popx<-samples[,c(8:7,2)] ; popfor<-aggregate(popx[,1:2], list(popx[,3]), mean) #long:lat, popn
bestlike<-which(likevalues[likevalues[,2]==j,1]==max(likevalues[likevalues[,2]==j,1]))
admix<-as.matrix(read.table(paste("./",prefix_sample_file,"_k",j,"/",prefix_sample_file,"_k",j,"_s",bestlike,".qopt",sep=""))) 
rownames(admix)<-seq(length=nrow(admix)) ;
admix<-admix[-c(urbsamples),]
admix<-aggregate(admix[,1:ncol(admix)], list(samples[,2]), sum); 
admix[,ncol(admix)+1]<-NA#admix[,ncol(admix)]
for(i in 1:length(admix[,1])){ admix[i,ncol(admix)]<-sum(admix[i,2:(ncol(admix)-1)]) }
popfor[,4:(ncol(popfor)+(ncol(admix)-1))]<-admix[,2:ncol(admix)]
colnames(popfor)[4:ncol(popfor)]<-c(LETTERS[1:(ncol(admix)-2)],"radius")
#popfor$radius <- (popfor$radius)
#plot 
p2 <- ggplot() + #data = rast.table, aes(x = x, y = y)
  #geom_tile(fill = rast.table$rgb) +
  #scale_alpha_discrete(range=c(1,0)) +
  geom_path(data=count.subset, aes(x = long, y = lat, group = group), color = 'black', size=0.5) +
  geom_path(data=river.subset, aes(x = long, y = lat, group = group), color = 'blue', alpha=0.7, size=0.3) +
  geom_path(data=road.subset, aes(x = long, y = lat, group = group), color = 'azure3', alpha=0.7, size=0.3) +
  geom_path(data=rail.subset, aes(x = long, y = lat, group = group), color = 'azure2', alpha=0.7, size=0.3) +
  geom_polygon(data=urban.subset, aes(x = long, y = lat, group = group), fill = "beige", colour = "beige") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('') + ylab('')
p2 <- p2 + geom_scatterpie(aes(x=V8, y=V7, r=0.25), data=popfor, cols=LETTERS[1:(ncol(admix)-2)], alpha=1, color=NA) #alpha is transparency for remove outline
#getting lineage informaiton onto dataframe
popfor <- cbind(popfor, lineage=samples[match(popfor[,1], samples[,2]),9])
p2 <- p2 + geom_text_repel(data=popfor, aes(x=V8, y=V7, label=Group.1, colour=lineage), min.segment.length = 0.1, size=1, segment.size=0.1, bg.color = "black", bg.r = .05) # geom_text_repel(data = popfor, aes(x = V8, y = V7, label = Group.1), size=2, colour="yellow") # fontface="bold",
p2 <- p2 + theme_bw() +
  theme(legend.position = "bottom",plot.margin = unit(c(0, 0.5, 0, 0), 'lines'),
        panel.border = element_rect(colour = 'black', fill = 'transparent')) + ggtitle(paste0("Without Urban; k=", j))
#p2 + geom_scatterpie_legend(c(.015,.030,.060), x=49, y=-16.6,n=3)
print(p2) #1*radius
}

  
dev.off()

##########################################
####### Plotting all samples normalised He ########## 
##########################################
pdf(paste("./",prefix_sample_file,"_geo_rural_He.pdf",sep=""), 10, 6) 
#par(mai = c(2, 2, 2, 5))
  samples<-read.table(popmap,header=F)
  #for X only 
   #samples<-samples[-which(samples$V1 == "Y199_1"),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
   urbsamples<-which(samples$V10 %in% c("Urban_Market", "Seizure")) # %in% c("Yaounde", "Abgab", "Tongue", "Benin", "DCM", "Cote_dIvoire", "DDM", "Ghana", "SeizureP", "SeizureB", "Central_African_Republic"))
   samples<-samples[-c(urbsamples),] ; samples<-droplevels(samples) ; rownames(samples) <- NULL
  popx<-samples[,c(8:7,2)] ; popfor<-aggregate(popx[,1:2], list(popx[,3]), mean) #long:lat, popn

p2 <- ggplot() + #data = rast.table, aes(x = x, y = y)
  #geom_tile(fill = rast.table$rgb) +
  #scale_alpha_discrete(range=c(1,0)) +
  geom_path(data=count.subset, aes(x = long, y = lat, group = group), color = 'black', size=0.5) +
  geom_path(data=river.subset, aes(x = long, y = lat, group = group), color = 'blue', alpha=0.7, size=0.3) +
  geom_path(data=road.subset, aes(x = long, y = lat, group = group), color = 'azure3', alpha=0.7, size=0.3) +
  geom_path(data=rail.subset, aes(x = long, y = lat, group = group), color = 'azure2', alpha=0.7, size=0.3) +
  geom_polygon(data=urban.subset, aes(x = long, y = lat, group = group), fill = "beige", colour = "beige") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('') + ylab('')
p2 <- p2 + geom_point(data=samples, aes(x=V8, y=V7, color=V12), alpha=0.5, size=3) + geom_text_repel(data=popfor, aes(x=V8, y=V7, label=Group.1), min.segment.length = 0.1, size=1, segment.size=0.1, bg.color = "black", bg.r = .05) + theme_bw() +
  theme(legend.position = "bottom",plot.margin = unit(c(0, 0.5, 0, 0), 'lines'),
        panel.border = element_rect(colour = 'black', fill = 'transparent'))  +  scale_color_continuous(type = "viridis")
        # for border around circles +  geom_point(data=samples, aes(x=V8, y=V7), shape = 1,size = 3,colour = "black") 
print(p2)
dev.off()




##########################################
####### just plot map ########## 
##########################################
rast.new = rast.table[seq(1, nrow(rast.table), 5), ] #this makes it ugly but far more loadable by only keeping every 5th line

pdf(paste("./",prefix_sample_file,"map.pdf",sep=""), 10, 6) 
p2 <- ggplot(data = rast.new, aes(x = x, y = y)) + 
  geom_tile(fill = rast.new$rgb) +
  scale_alpha_discrete(range=c(1,0)) +
  geom_path(data=count.subset, aes(x = long, y = lat, group = group), color = 'black', size=0.5) +
  geom_path(data=river.subset, aes(x = long, y = lat, group = group), color = 'deepskyblue1', alpha=0.7, size=0.3) +
  geom_path(data=road.subset, aes(x = long, y = lat, group = group), color = 'azure4', alpha=0.7, size=0.3) +
  geom_path(data=rail.subset, aes(x = long, y = lat, group = group), color = 'azure3', alpha=0.7, size=0.3) +
  geom_polygon(data=urban.subset, aes(x = long, y = lat, group = group), fill = "darkgray", colour = "darkgray") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('') + ylab('')
print(p2)
  
dev.off()


q(save="no")



