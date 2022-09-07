#!/bin/sh
#!/usr/bin/env Rscript

#module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R
#Rscript $bin/plot_locator.R $loc_dir/ allsamples $res/locator/metadata.txt $res/locator/unknownsamples.txt $loc

#--sample_data $res/locator/metadata.txt --out $dir/alltest_ --map T --longlat T

library(ggplot2)
library(rgdal)
library(plyr)
library(dplyr)
library(raster)
library(ggnewscale)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)


loc_dir <- args[1]
pred_file <- args[2]
metadata <- args[3]
oldloc <- args[4]
loc <- args[5]

print(c(loc_dir, pred_file, metadata, oldloc, loc))
#for non-iteration
#loc_dir <- "/work/idumville/pangolins/RAD/results/locator.xchrom"
#pred_file <- "allsamples"
#metadata <- "/work/idumville/pangolins/RAD/results/locator/metadata.txt"
#oldloc <- "/work/idumville/pangolins/RAD/results/locator/unknownsamples.txt"
#loc <- "xchrom"

setwd(loc_dir)

#locator_dir="/work/idumville/pangolins/RAD/results/locator.xchrom"


predicted<-read.table(paste("./",pred_file,"_predlocs.txt",sep=""),header=T,sep=",")
predicted$t<-"predicted"
known<-read.table(metadata,header=T)
known<-na.omit(known)
known$t<-"training"
oldloc<-read.table(oldloc,header=T)
oldloc$t<-"collection"
traced<-rbind(known, oldloc)
traced<-rbind(traced, predicted)

map <- readRDS("/work/idumville/pangolins/RAD/data/geodata/bigmap.rdata")


pdf(paste("./",loc,"_locator_overallresults.pdf",sep=""), 10, 6)

p1 <- map + new_scale_fill() + scale_fill_brewer(palette = "Dark2") + scale_color_brewer(palette = "Dark2") + geom_point(data=traced, mapping=aes(x=y, y=x, fill=t), size=0.6, colour="black", pch=21) + geom_point(data=traced[which(traced$t %in% c("training", "collection")),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black",pch=21) + theme(legend.position = "bottom") + geom_line(data=traced, aes(x=y, y=x, group = sampleID, alpha=0.5, size=0.5),color="darkslategrey", size=0.1) + theme_bw() #map + xlim(lon) + ylim(lat) + 
plot(p1)
dev.off()

map <- readRDS("/work/idumville/pangolins/RAD/data/geodata/smallmap.rdata")

lon <- c(8.000, 14.2194) #-10.5318, 29.53505)
lat <- c(-0.271891, 7.25)

pdf(paste("./",loc,"_locator_cameroonresults.pdf",sep=""), 6, 6)

p1 <- map + new_scale_fill() + scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") + 
geom_point(data=traced, mapping=aes(x=y, y=x, fill=t), size=0.6, colour="black",pch=21) + 
geom_point(data=traced[which(traced$t %in% c("training", "collection")),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black", pch=21) + 
theme(legend.position = "bottom") + 
geom_line(data=traced, aes(x=y, y=x, group = sampleID, alpha=0.5, size=0.5),color="darkslategrey", size=0.1) + 
theme_bw() + xlim(lon) + ylim(lat)
plot(p1)
dev.off()

 loc_dir <- "/work/idumville/pangolins/RAD/results/locator"
pred_file <- "/work/idumville/pangolins/RAD/results/locator/indvPredictedLong.txt"
metadata <- "/work/idumville/pangolins/RAD/results/locator/metadata.txt"
oldloc <- "/work/idumville/pangolins/RAD/results/locator/unknownsamples.txt"
setwd(loc_dir)

known<-read.table(metadata,header=T)
known<-na.omit(known)
known$t<-"training"
oldloc<-read.table(oldloc,header=T)
oldloc$t<-"collection"
traced<-rbind(known, oldloc)

predicted<-read.table(pred_file, header=T, col.names=c("sampleID", "y", "x", "t"))
predicted <- predicted[which(predicted$t == "kd"),]
traced<-rbind(traced, predicted) 
library(stringr)
traced[c('Indv', 'Rep')] <- str_split_fixed(traced$sampleID, "_", 2)



pdf(paste(loc_dir, "/allindv.pdf",sep=""), 6, 6)
p1 <- map + new_scale_fill() + scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") + 
geom_point(data=traced, mapping=aes(x=y, y=x, fill=t), size=0.6, colour="black",pch=21) + 
geom_point(data=traced[which(traced$t %in% c("training", "collection")),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black", pch=21) + 
theme(legend.position = "bottom") + 
geom_line(data=traced, aes(x=y, y=x, group = sampleID), alpha =0.5,color="darkslategrey", size=0.1) + 
theme_bw() + xlim(lon) + ylim(lat) + scale_fill_discrete(name = NULL, labels = c("Collection", "Peak Kernal Density", "Training"))
plot(p1)
dev.off()

predicted[c('Indv', 'Rep')] <- str_split_fixed(predicted$sampleID, "_", 2)
repeatedsamps <- predicted[which(predicted$Rep == 2), 5]

reppred <- predicted[which(predicted$Indv %in% repeatedsamps),] 
reppredkd <- reppred[which(reppred$t == "kd"),] 

pdf(paste(loc_dir, "/repeatedindvkd.pdf",sep=""), 6, 6)
p1 <- map + new_scale_fill() + scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") + 
geom_point(data=reppredkd, mapping=aes(x=y, y=x, fill=Rep), size=1, colour="black",pch=21) + 
#geom_point(data=traced[which(traced$t %in% c("training", "collection")),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black", pch=21) + 
theme(legend.position = "bottom") + 
geom_line(data=reppredkd, aes(x=y, y=x, group = Indv), alpha =0.5,color="darkslategrey", size=0.1) + 
theme_bw() + xlim(lon) + ylim(lat) #+ scale_fill_discrete(name = NULL, labels = c("Collection", "Peak Kernal Density", "Training"))
plot(p1)
dev.off()

##Final Repeat Plot
png(paste("./minmac3_locator_repeats.png",sep=""),width=4, height=4.7,res = 600,units = "in")
p1 <- map + new_scale_fill() +  
scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") + 
geom_point(data=reppredkd, mapping=aes(x=y, y=x, fill=Rep), size=1.2, colour="black",pch=21) + 
geom_line(data=reppredkd, aes(x=y, y=x, group = Indv),color="darkslategrey", size=0.1) + theme_bw() + theme(legend.title=element_blank(), legend.position = c(0.92,0.90), legend.key=element_blank(), legend.background=element_blank()) + labs(tag="A") +
   theme(plot.margin = unit(c(0.1,-0.5,-0.3,-1), "cm"), plot.tag.position = c(0.04, 0.98))
plot(p1)
dev.off()



### for 00 and 90%

 loc_dir <- "/work/idumville/pangolins/RAD/results/locator"
pred_file <- "/work/idumville/pangolins/RAD/results/locator/seeding.concat.allWCA/seedinglocator90/seeding90error_centroids.txt"
metadata <- "/work/idumville/pangolins/RAD/results/locator/metadata.txt"
oldloc <- "/work/idumville/pangolins/RAD/results/locator/unknownsamples.txt"

known<-read.table(metadata,header=T)
known<-na.omit(known)
known$t<-"training"
oldloc<-read.table(oldloc,header=T)
oldloc$t<-"collection"
traced<-rbind(known, oldloc)

predicted<-read.table(pred_file, header=T, col.names=c("sampleID", "y", "x", "gdy", "gdx"))
predicted <- predicted[,-c(4,5)]
predicted$t <- "kd"
traced<-rbind(traced, predicted) 
library(stringr)
traced[c('Indv', 'Rep')] <- str_split_fixed(traced$sampleID, "_", 2)


map <- readRDS("/work/idumville/pangolins/RAD/data/geodata/smallmap.rdata")

lon <- c(8.000, 14.2194) #-10.5318, 29.53505)
lat <- c(-0.271891, 7.25)

p1 <- map + new_scale_fill() +# scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") + 
geom_point(data=traced, mapping=aes(x=y, y=x, fill=t), size=0.6, colour="black",pch=21) + 
geom_point(data=traced[which(traced$t %in% c("training", "collection")),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black", pch=21) + theme_bw()  +
  theme(legend.title=element_blank(), legend.position = c(0.86,0.94), legend.key=element_blank(), legend.background=element_blank()) + labs(tag="A") +
   theme(plot.tag.position = c(0.04, 0.98)) + labs(tag="B") +
geom_line(data=traced, aes(x=y, y=x, group = sampleID), alpha =0.5,color="darkslategrey", size=0.1) + 
 scale_fill_manual(values = c( "#1B9E77", "#D95F02" ,"#7570B3"), labels = c("Collection", "Prediction", "Training")) +
coord_cartesian(
  xlim = lon,
  ylim = lat)


loc_dir <- "/work/idumville/pangolins/RAD/results/locator"
pred_file <- "/work/idumville/pangolins/RAD/results/locator/seeding.concat.allWCA/seedinglocator00/seeding00error_centroids.txt"
metadata <- "/work/idumville/pangolins/RAD/results/locator/metadata.txt"
oldloc <- "/work/idumville/pangolins/RAD/results/locator/unknownsamples.txt"

known<-read.table(metadata,header=T)
known<-na.omit(known)
known$t<-"training"
oldloc<-read.table(oldloc,header=T)
oldloc$t<-"collection"
traced<-rbind(known, oldloc)

predicted<-read.table(pred_file, header=T, col.names=c("sampleID", "y", "x", "gdy", "gdx"))
predicted <- predicted[,-c(4,5)]
predicted$t <- "kd"
traced<-rbind(traced, predicted) 
library(stringr)
traced[c('Indv', 'Rep')] <- str_split_fixed(traced$sampleID, "_", 2)

p2 <- map + new_scale_fill() +# scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") +  
 scale_fill_manual(values = c( "#1B9E77", "#D95F02" ,"#7570B3"), labels = c("Collection", "Prediction", "Training")) +
geom_point(data=traced, mapping=aes(x=y, y=x, fill=t), size=0.6, colour="black",pch=21) + 
geom_point(data=traced[which(traced$t %in% c("training", "collection")),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black", pch=21) +  theme_bw() +
  theme(legend.position = "none")  + labs(tag="A") +
   theme(plot.tag.position = c(0.04, 0.98)) +
        labs(tag="A") +
geom_line(data=traced, aes(x=y, y=x, group = sampleID), alpha =0.5,color="darkslategrey", size=0.1) + 
 coord_cartesian( xlim = lon,  ylim = lat)
 
require(gridExtra)

png(paste(loc_dir, "/Supp_S18_00and90allindv.png",sep=""), width = 8, height = 5, res=600, units="in")
grid.arrange(p2, p1, ncol=2)
dev.off()





## repeats
#90
pred_file <- "/work/idumville/pangolins/RAD/results/locator/seeding.concat.allWCA/seedinglocator90/seeding90error_centroids.txt"
oldloc <- "/work/idumville/pangolins/RAD/results/locator/unknownsamples.txt"
known<-read.table(metadata,header=T)
known<-na.omit(known)
known$t<-"training"
oldloc<-read.table(oldloc,header=T)
oldloc$t<-"collection"
traced<-rbind(known, oldloc)
predicted<-read.table(pred_file, header=T, col.names=c("sampleID", "y", "x", "gdy", "gdx"))
predicted <- predicted[,-c(4,5)]
predicted$t <- "kd"
traced<-rbind(traced, predicted) 
library(stringr)
traced[c('Indv', 'Rep')] <- str_split_fixed(traced$sampleID, "_", 2)
predicted[c('Indv', 'Rep')] <- str_split_fixed(predicted$sampleID, "_", 2)
repeatedsamps <- predicted[which(predicted$Rep == 2), 5]
reppred <- predicted[which(predicted$Indv %in% repeatedsamps),] 
reppredkd <- reppred[which(reppred$t == "kd"),] 

rep90 <- map + new_scale_fill() + scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") + 
geom_point(data=reppredkd, mapping=aes(x=y, y=x, fill=Rep), size=1, colour="black",pch=21) + 
#geom_point(data=traced[which(traced$t %in% c("training", "collection")),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black", pch=21) + 
geom_line(data=reppredkd, aes(x=y, y=x, group = Indv), color="darkslategrey", size=0.15) + 
theme_bw() +    theme(    legend.position="none",
        axis.title = element_text(size=6),plot.tag.position = c(0.04, 0.98)) + labs(tag="B") +
        coord_cartesian(  xlim = lon,  ylim = lat)

#00
pred_file <- "/work/idumville/pangolins/RAD/results/locator/seeding.concat.allWCA/seedinglocator00/seeding00error_centroids.txt"
oldloc <- "/work/idumville/pangolins/RAD/results/locator/unknownsamples.txt"
known<-read.table(metadata,header=T)
known<-na.omit(known)
known$t<-"training"
oldloc<-read.table(oldloc,header=T)
oldloc$t<-"collection"
traced<-rbind(known, oldloc)
predicted<-read.table(pred_file, header=T, col.names=c("sampleID", "y", "x", "gdy", "gdx"))
predicted <- predicted[,-c(4,5)]
predicted$t <- "kd"
traced<-rbind(traced, predicted) 
library(stringr)
traced[c('Indv', 'Rep')] <- str_split_fixed(traced$sampleID, "_", 2)
predicted[c('Indv', 'Rep')] <- str_split_fixed(predicted$sampleID, "_", 2)
repeatedsamps <- predicted[which(predicted$Rep == 2), 5]
reppred <- predicted[which(predicted$Indv %in% repeatedsamps),] 
reppredkd <- reppred[which(reppred$t == "kd"),] 
rep00 <- map + new_scale_fill() + scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") + 
geom_point(data=reppredkd, mapping=aes(x=y, y=x, fill=Rep), size=1, colour="black",pch=21) + 
geom_line(data=reppredkd, aes(x=y, y=x, group = Indv), color="darkslategrey", size=0.15) +   
theme(legend.position="none", plot.tag.position = c(0.04, 0.98)) +
        labs(tag="A") +
 coord_cartesian(  xlim = lon,  ylim = lat)
 
  loc_dir <- "/work/idumville/pangolins/RAD/results/locator"
pred_file <- "/work/idumville/pangolins/RAD/results/locator/indvPredictedLong.txt"
metadata <- "/work/idumville/pangolins/RAD/results/locator/metadata.txt"
oldloc <- "/work/idumville/pangolins/RAD/results/locator/unknownsamples.txt"
setwd(loc_dir)

known<-read.table(metadata,header=T)
known<-na.omit(known)
known$t<-"training"
oldloc<-read.table(oldloc,header=T)
oldloc$t<-"collection"
traced<-rbind(known, oldloc)

predicted<-read.table(pred_file, header=T, col.names=c("sampleID", "y", "x", "t"))
predicted <- predicted[which(predicted$t == "kd"),]
traced<-rbind(traced, predicted) 
library(stringr)
traced[c('Indv', 'Rep')] <- str_split_fixed(traced$sampleID, "_", 2)

predicted[c('Indv', 'Rep')] <- str_split_fixed(predicted$sampleID, "_", 2)
repeatedsamps <- predicted[which(predicted$Rep == 2), 5]

reppred <- predicted[which(predicted$Indv %in% repeatedsamps),] 
reppredkd <- reppred[which(reppred$t == "kd"),] 

repindv <- map + new_scale_fill() +  
scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") + 
geom_point(data=reppredkd, mapping=aes(x=y, y=x, fill=Rep), size=1.2, colour="black",pch=21) + 
geom_line(data=reppredkd, aes(x=y, y=x, group = Indv),color="darkslategrey", size=0.15) + theme_bw() + theme(legend.title=element_blank(),  legend.position = c(0.90,0.90), legend.key=element_blank(), legend.background=element_blank()) + labs(tag="C") +
   theme(plot.tag.position = c(0.04, 0.98)) +
 coord_cartesian(  xlim = lon,  ylim = lat)


map <- readRDS("/work/idumville/pangolins/RAD/data/geodata/smallmap_noplane.rdata")


png(paste(loc_dir, "/Supp_S19_allrepeats.png",sep=""), width = 12, height = 5, res=600, units="in")
grid.arrange(rep00, rep90, repindv, ncol=3)
dev.off()


############################################
######## Final AL plot #####################
##############################################
# %in% c("Yaounde_Nkolndongo_Market", "Douala_Dakat_Market", "SeizureR", "SeizureP", "Douala_Central_Market", "SeizureB", "Banyo"),]


infofile="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/sampledataautosomes.txt"
samples<-read.table(infofile,header=F,sep="\t")
samplestorem <- samples[which(samples$V2 %in% c("Banyo")), "V1"]

loc_dir <- "/work/idumville/pangolins/RAD/results/locatorAL"
pred_file <- "/work/idumville/pangolins/RAD/results/locatorAL/indvlocmm2predlocs.txt"
metadata <- "/work/idumville/pangolins/RAD/results/locatorAL/metadata.txt"
oldloc <- "/work/idumville/pangolins/RAD/results/locatorAL/unknownsamples.txt"

loc="ALindvmm2"

setwd(loc_dir)


known<-read.table(metadata,header=T)
known<-na.omit(known)
known$t<-"Training"
oldloc<-read.table(oldloc,header=T)
oldloc$t<-"Collection"
traced<-rbind(known, oldloc)

predicted<-read.table(pred_file, header=T, col.names=c("sampleID", "y", "x", "gdy", "gdx"))
predicted <- predicted[,-c(4,5)]
predicted$t <- "Predicted"
traced<-rbind(traced, predicted) 
library(stringr)
traced[c('Indv', 'Rep')] <- str_split_fixed(traced$sampleID, "_", 2)

map <- readRDS("/work/idumville/pangolins/RAD/data/geodata/bigmap.rdata")

traced$y <- as.numeric(traced$y)
traced$x <- as.numeric(traced$x)

traced <- traced[-which(traced$sampleID == "DlaB89_1"),]
traced <- traced[-which(traced$sampleID %in% samplestorem),] #removing Banyo

traced[which(traced$y == "8.5" & traced$t == "Collection"), "y"] <- 4.000 #moving seizures

#pdf(paste("./Figure5A",loc,"_locator_overallresults.pdf",sep=""), 10, 4)
png(paste("./Figure5A",loc,"_locator_overallresults.png",sep=""),width=7, height=3,res = 600,units = "in")
p1 <- map + new_scale_fill() +  
scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") + 
geom_point(data=traced[which(traced$t %in% c("Training", "Collection")),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black",pch=21) +
 geom_point(data=traced[which(traced$t  == "Predicted"),], mapping=aes(x=y, y=x, fill=t), size=1.2, colour="black", pch=21)   + 
 geom_line(data=traced, aes(x=y, y=x, group = sampleID), alpha=0.5,color="darkslategrey", size=0.1) + theme_bw() + theme(legend.title=element_blank(), legend.position = c(0.07,0.2), legend.key=element_blank(), legend.background=element_blank()) + labs(tag="A") +
   theme(plot.margin = unit(c(0,0.1,-0.7,-0.4), "cm"), plot.tag.position = c(0.03, 0.98))
plot(p1)
dev.off()


#just plot seizures
traced <- traced[which(traced$sampleID %in% c("LL180926-003_1", "LL180926-022-A_1" , "LL180926-029-A_1",
"T-1654_1", "T-2139_1" , "T-2145_1" ,"T-2163_1", "T-3349_1", "T-3350_1")),]

traced <- subset(traced, y < 14.2194 & y > 8)
traced <- subset(traced, x < 7.25 & x > -0.271891)


png(paste("./",loc,"_locator_seizure.png",sep=""),width=7, height=3,res = 600,units = "in")
p1 <- map + new_scale_fill() +  
scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") + 
geom_point(data=traced[which(traced$t %in% c("Training", "Collection")),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black",pch=21) +
 geom_point(data=traced[which(traced$t  == "Predicted"),], mapping=aes(x=y, y=x, fill=t), size=1.2, colour="black", pch=21)   + 
 geom_line(data=traced, aes(x=y, y=x, group = sampleID), alpha=0.5,color="darkslategrey", size=0.1) + theme_bw() + theme(legend.title=element_blank(), legend.position = c(0.07,0.4), legend.key=element_blank(), legend.background=element_blank()) #map + xlim(lon) + ylim(lat) + 
plot(p1)
dev.off()

############################################
######## Final WCA plot #####################
##############################################


loc_dir <- "/work/idumville/pangolins/RAD/results/locator"
pred_file <- "/work/idumville/pangolins/RAD/results/locator/indvPredictedLong.txt"
metadata <- "/work/idumville/pangolins/RAD/results/locator/metadata.txt"
oldloc <- "/work/idumville/pangolins/RAD/results/locator/unknownsamples.txt"

loc="WCAindvmm3"

setwd(loc_dir)

known<-read.table(metadata,header=T)
known<-na.omit(known)
known$t<-"Training"
oldloc<-read.table(oldloc,header=T)
oldloc$t<-"Collection"
traced<-rbind(known, oldloc)

predicted<-read.table(pred_file, header=T, col.names=c("sampleID", "y", "x", "kd"))
predicted <- predicted[which(predicted$kd=="kd"),]
predicted$t <- "Predicted"
predicted <- predicted[,c(1,2,3,5)]
traced<-rbind(traced, predicted) 
library(stringr)
#traced[c('Indv', 'Rep')] <- str_split_fixed(traced$sampleID, "_", 2)

map <- readRDS("/work/idumville/pangolins/RAD/data/geodata/smallmap.rdata")

traced$y <- as.numeric(traced$y)
traced$x <- as.numeric(traced$x)

#lon <- c(8.000, 14.2194)
#lat <- c(-0.271891, 7.25) 

traced <- subset(traced, y < 14.2194 & y > 8)
traced <- subset(traced, x < 7.25 & x > -0.271891)


#adding size of reference markets

meta <- read.table("/work/idumville/pangolins/RAD/results/allmarkets.txt", header=T)
marks <- meta[,c(4,5)]
ag <- aggregate(.~Locality, meta, FUN=head, 1)
count <- meta %>% count(Locality)
ag$count <- count$n
ag <- distinct(merge(ag, marks, by="Locality", all.y=FALSE))

ag$Latiude <- as.numeric(ag$Latitude)
ag$Longitude <- as.numeric(ag$Longitude)

ag <- ag[which(ag$Market_Type.y == "Rural_Market"),]
ag$Market_Type.y <- "Training" 
ag <- subset(ag, Longitude < 14.2194 & Longitude > 8)
#ag <- subset(ag, Latitude < 7.25 & Latitude > -0.271891)


#pdf(paste("./",loc,"_locator_overallresults.pdf",sep=""), 10, 4)
png(paste("./",loc,"_locator_overallresultstemp.png",sep=""),width=4, height=4.7,res = 600,units = "in")
p1 <- map + new_scale_fill() +  scale_fill_brewer(palette = "Dark2") + scale_color_brewer(palette = "Dark2") +
 geom_point(data=traced[which(traced$t %in% c("Training", "Collection")),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black",pch=21) + 
geom_point(data=traced[which(traced$t  == "Predicted"),], mapping=aes(x=y, y=x, fill=t), size=1.2, colour="black", pch=21)   + 
geom_line(data=traced, aes(x=y, y=x, group = sampleID), alpha=0.5,color="darkslategrey", size=0.1) + theme_bw() + guides(fill = guide_legend(label.position = "left")) + theme(legend.title=element_blank(), legend.position = c(0.88,0.93), legend.key=element_blank(), legend.background=element_blank()) 
plot(p1)
dev.off()


png(paste("./Figure6A",loc,"_locator_overallresults2.png",sep=""),width=4, height=4.7,res = 600,units = "in")
p1 <- map + new_scale_fill() +  scale_fill_brewer(palette = "Dark2") + scale_color_brewer(palette = "Dark2") +
 geom_point(data=traced[which(traced$t == "Collection"),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black",pch=21) + 
geom_point(data=ag, mapping=aes(x=Longitude, y=Latitude, fill=Market_Type.y, size=count), colour="black", pch=21)   + 
geom_point(data=traced[which(traced$t  == "Predicted"),], mapping=aes(x=y, y=x, fill=t), size=1.2, colour="black", pch=21)   + 
geom_line(data=traced, aes(x=y, y=x, group = sampleID), alpha=0.5,color="darkslategrey", size=0.1) + theme_bw() + theme(legend.title=element_blank(), legend.position = c(0.3,0.95), legend.spacing = unit(0,'pt'), legend.margin = margin(t=0,b=0,unit='pt'), legend.key=element_blank(), legend.background=element_blank(), legend.direction="horizontal") +
  theme(legend.key.size = unit(0.1, 'in'), #change legend key size
        legend.key.height = unit(0.1, 'in'), #change legend key height
        legend.key.width = unit(0.1, 'in'), #change legend key width
        legend.text = element_text(size=8)) + labs(tag="A") +
  theme(plot.margin = unit(c(0.1,-0.5,-0.3,-1), "cm"), plot.tag.position = c(0.05, 0.98))
plot(p1)
dev.off()


traced <- traced[which(traced$sampleID %in% c("LL180926-003_1", "LL180926-022-A_1" , "LL180926-029-A_1",
"T-1654_1", "T-2139_1" , "T-2145_1" ,"T-2163_1", "T-3349_1", "T-3350_1")),]
png(paste("./",loc,"_locator_WCAseizure.png",sep=""),width=4, height=4.7,res = 600,units = "in")
p1 <- map + new_scale_fill() +  scale_fill_brewer(palette = "Dark2") + scale_color_brewer(palette = "Dark2") +
 geom_point(data=traced[which(traced$t == "Collection"),], mapping=aes(x=y, y=x, fill=t), size=2, colour="black",pch=21) + 
geom_point(data=ag, mapping=aes(x=Longitude, y=Latitude, fill=Market_Type.y, size=count), colour="black", pch=21)   + 
geom_point(data=traced[which(traced$t  == "Predicted"),], mapping=aes(x=y, y=x, fill=t), size=1.2, colour="black", pch=21)   + 
geom_line(data=traced, aes(x=y, y=x, group = sampleID), alpha=0.5,color="darkslategrey", size=0.1) + theme_bw() + theme(legend.title=element_blank(), legend.position = c(0.35,0.05), legend.spacing = unit(0,'pt'), legend.margin = margin(t=0,b=0,unit='pt'), legend.key=element_blank(), legend.background=element_blank(), legend.direction="horizontal") +
  theme(legend.key.size = unit(0.1, 'in'), #change legend key size
        legend.key.height = unit(0.1, 'in'), #change legend key height
        legend.key.width = unit(0.1, 'in'), #change legend key width
        legend.text = element_text(size=8))
plot(p1)
dev.off()

