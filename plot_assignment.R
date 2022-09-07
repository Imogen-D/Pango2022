

#module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R


library(ggplot2)
library(rgdal)
library(plyr)
library(dplyr)
library(raster)
library(ggnewscale)
library(RColorBrewer)

dir <- "/work/idumville/pangolins/RAD/results/assignment"
setwd(dir)


BONEpredicted<-read.delim("BONEcoords.txt", header=T) #made in PlottingAssignment.R (local)
BONEpredicted$t <- ifelse(BONEpredicted$t == "Assigned", "Predicted", ifelse(BONEpredicted$t == "Origin", "Collection", "Reference"))


map <- readRDS("/work/idumville/pangolins/RAD/data/geodata/bigmap.rdata")


BONEpredicted[which(BONEpredicted$Locality %in% c("SeizureB", "SeizureP", "SeizureR")), "Longitude"] <- 4.000

p1 <- map + new_scale_fill() +  scale_fill_brewer(palette = "Dark2") + scale_color_brewer(palette = "Dark2") +
geom_point(data=BONEpredicted[which(BONEpredicted$t == "Reference"),], mapping=aes(x=Longitude, y=Latitude, fill=t), size=2, colour="black",pch=21) +
geom_point(data=BONEpredicted[which(BONEpredicted$t == "Predicted"),], mapping=aes(x=Longitude, y=Latitude, fill=t), size=2, colour="black",pch=21) +  
geom_point(data=BONEpredicted[which(BONEpredicted$t == "Collection"),], mapping=aes(x=Longitude, y=Latitude, fill=t), size=2, colour="black",pch=21) + 
geom_line(data=BONEpredicted[which(BONEpredicted$t != "Reference"),], aes(x=Longitude, y=Latitude, group = new_name), alpha=0.5, size=0.16,color="darkslategrey") + theme_bw() + theme(legend.title=element_blank(), legend.position = c(0.07,0.2), legend.key=element_blank(), legend.background=element_blank()) 

png("Supp_FigureS15_ALBONEresults.png",width=7, height=3,res = 600,units = "in")
plot(p1)
dev.off()


BONEpredicted<-read.delim("BONEcoords.txt", header=T)
BONEpredicted$t <- ifelse(BONEpredicted$t == "Assigned", "Predicted", ifelse(BONEpredicted$t == "Origin", "Collection", "Reference"))

BONEpredicted <- BONEpredicted[-which(BONEpredicted$Locality == "Banyo"),]

smolmap <- readRDS("/work/idumville/pangolins/RAD/data/geodata/smallmap.rdata")

BONEpredicted <- subset(BONEpredicted, Longitude < 14.2194 & Longitude > 8)
BONEpredicted <- subset(BONEpredicted, Latitude < 7.25 & Latitude > -0.271891)

BONEassign <- BONEpredicted[which(BONEpredicted$t == "Predicted"),]
#just getting top and then removing duplicates - not ideal way to do it as just removes the seocnd one
BONEassign <- BONEassign %>% group_by(new_name) %>% top_n(1, value)
BONEassign <- BONEassign[-c(which(duplicated(BONEassign$new_name))),]
BONEassign <- merge(data.frame(BONEassign), as.data.frame(BONEassign) %>% count(Locality))

BONEpredicted$n <- "NA"

BONEref <- BONEpredicted[which(BONEpredicted$t == "Reference"),]
BONEcollect <- BONEpredicted[which(BONEpredicted$t == "Collection"),]
#just getting top and then removing duplicates - not ideal way to do it as just removes the seocnd one
BONEcollect <- BONEcollect %>% group_by(new_name) %>% top_n(1, value)
BONEcollect <- BONEcollect[-c(which(duplicated(BONEcollect$new_name))),]


BONEpredicted <- rbind(data.frame(BONEcollect), data.frame(BONEassign))
BONEpredicted <- rbind(data.frame(BONEpredicted), data.frame(BONEref))

p2 <- smolmap + new_scale_fill() +  scale_fill_brewer(palette = "Dark2") + scale_color_brewer(palette = "Dark2") +
geom_point(data=BONEpredicted[which(BONEpredicted$t == "Reference"),], mapping=aes(x=Longitude, y=Latitude, fill=t), size=2, colour="black",pch=21) + 
geom_point(data=BONEassign, mapping=aes(x=Longitude, y=Latitude, fill=t, size=n), colour="black",pch=21) + 
geom_point(data=BONEpredicted[which(BONEpredicted$t == "Collection"),], mapping=aes(x=Longitude, y=Latitude, fill=t), size=2, colour="black",pch=21) + 
geom_line(data=BONEpredicted[which(BONEpredicted$t != "Reference"),], aes(x=Longitude, y=Latitude, group = new_name), alpha=0.5, size=0.16,color="darkslategrey") + theme_bw() + theme(legend.title=element_blank(), legend.position = c(0.9,0.85), legend.key=element_blank(), legend.background=element_blank())  + labs(tag="A") + theme(plot.margin = unit(c(0.1,-0.5,-0.3,-1), "cm"), plot.tag.position = c(0.04, 0.98))

png("Fig4A_WCABONEresults.png",width=5, height=6,res = 400,units = "in")
plot(p2)
dev.off()

BONEpredicted[which(BONEpredicted$t != "Reference"),] %>% group_by(new_name)


###Supp_FigS14 comparing BONE

BONEpredictedWTA<-read.delim("BONEcoords_wWTA.txt", header=T)
BONEpredictedWTA$t <- ifelse(BONEpredictedWTA$t == "Assigned", "Predicted", ifelse(BONEpredictedWTA$t == "Origin", "Collection", ifelse(BONEpredictedWTA$t == "WTA", "WTA", "Reference")))

smolmap <- readRDS("/work/idumville/pangolins/RAD/data/geodata/smallmap.rdata")

BONEpredictedWTA <- subset(BONEpredictedWTA, Longitude < 14.2194 & Longitude > 8)
BONEpredictedWTA <- subset(BONEpredictedWTA, Latitude < 7.25 & Latitude > -0.271891)

BONEassignWTA <- BONEpredictedWTA[which(BONEpredictedWTA$t %in% c("Predicted", "WTA")),]

#BONEassignWTA <- merge(BONEassignWTA, BONEassignWTA %>% count(Locality))

BONEassignWTA <- BONEassignWTA[-which(BONEassignWTA$Locality == "Banyo"),]
BONEpredictedWTA <- BONEpredictedWTA[-which(BONEpredictedWTA$Locality == "Banyo"),]

WTAassign <- merge(BONEassignWTA[which(BONEassignWTA$t =="WTA"),], BONEassignWTA[which(BONEassignWTA$t =="WTA"),] %>% count(Locality))
WTAassign <- WTAassign[match(unique(WTAassign$new_name), WTAassign$new_name),]

SPassign <- merge(BONEassignWTA[which(BONEassignWTA$t =="Predicted"),], BONEassignWTA[which(BONEassignWTA$t =="Predicted"),] %>% count(Locality))
SPassign <- SPassign[match(unique(SPassign$new_name), SPassign$new_name),]

length(which(paste(SPassign$new_name, SPassign$Locality) %in% paste(WTAassign$new_name, WTAassign$Locality)))
length(which(paste(WTAassign$new_name, WTAassign$Locality) %in% paste(SPassign$new_name, SPassign$Locality)))

BONEassignWTA <- rbind(WTAassign, SPassign)


png("Supp_FigS14_WCABONEresults.png",width=4, height=4,res = 400,units = "in")
#plot(p2)
ggplot(BONEassignWTA, aes(x=Locality, fill=t)) +  geom_bar(stat="count", position=position_dodge(), colour="black", size=0.2) + theme_bw() +
theme(legend.title = element_blank(),
        legend.position=c(x=0.1, y=0.85),
        axis.title = element_text(size=4),
        axis.text = element_text(size=3),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'), #change legend key width
        legend.text = element_text(size=3),
        legend.background=element_blank())       + ylab("Number of samples assigned") + xlab("") +
        scale_fill_brewer(palette="Set2", labels = c("SolPath",  "WTA"))
dev.off()

## getting locator info

loc_dir <- "/work/idumville/pangolins/RAD/results/locator"
pred_file <- "/work/idumville/pangolins/RAD/results/locator/indvPredictedLong.txt"
metadata <- "/work/idumville/pangolins/RAD/results/locator/metadata.txt"
oldloc <- "/work/idumville/pangolins/RAD/results/locator/unknownsamples.txt"

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

traced$y <- as.numeric(traced$y)
traced$x <- as.numeric(traced$x)

traced <- subset(traced, y < 14.2194 & y > 8)
traced <- subset(traced, x < 7.25 & x > -0.271891)

# getting NGSadmix info
NGSadmix <- read.table("/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads/assignmentofurban.txt", header=T)
#gab1 medoueu 2 Yaounde3 north4 bioko5 south6
NGSadmix$x <- 3.268185
NGSadmix$y <- 12.480171
NGSadmix$x[which(NGSadmix$assign == 5)] <- 3.517968 
NGSadmix$y[which(NGSadmix$assign == 5)]  <- 8.715969
NGSadmix$x[which(NGSadmix$assign == 1)] <- 1.580706 
NGSadmix$y[which(NGSadmix$assign == 1)]  <- 11.632576
NGSadmix$x[which(NGSadmix$assign == 2)] <- 0.994383
NGSadmix$y[which(NGSadmix$assign == 2)]  <- 10.7805982
NGSadmix$x[which(NGSadmix$assign == 3)] <- 3.8302938
NGSadmix$y[which(NGSadmix$assign == 3)]  <- 11.3703656
NGSadmix$x[which(NGSadmix$assign == 4)] <- 5.584278
NGSadmix$y[which(NGSadmix$assign == 4)]  <- 10.934584

NGSadmix <- merge(NGSadmix, NGSadmix %>% count(assign))
NGSadmix$t <- "Admix"

BONEassign$t <- "BONE"

loctraced <- traced[which(traced$t  == "Predicted"),]
loctraced$t <- "Locator"
loctraced
BONEassign2 <- BONEassign[,-c(1,3,7)]
names(BONEassign2) <- names(loctraced)
mix <- rbind(loctraced, BONEassign2)


smolmap <- readRDS("/work/idumville/pangolins/RAD/data/geodata/smallmap_noplane.rdata")

### plotting the BONE assignment, locator and NGSadmix on same?
p2 <- smolmap + new_scale_fill() +  scale_fill_brewer(palette = "Dark2") + scale_color_brewer(palette = "Dark2") + 
#geom_point(data=NGSadmix, mapping=aes(x=y, y=x, fill=t, size=n), colour="black",pch=21) +
geom_point(data=BONEassign, mapping=aes(x=Longitude, y=Latitude, fill=t, size=n), colour="black",pch=21)  +
geom_point(data=loctraced, mapping=aes(x=y, y=x, fill=t), size=1.2, colour="black", pch=21)  +
geom_line(data=mix, aes(x=y, y=x, group = sampleID), alpha=0.5, size=0.16,color="darkslategrey") +
labs(tag="A") + theme(plot.margin = unit(c(0.1,-0.5,-0.3,-1), "cm"), plot.tag.position = c(0.04, 0.98), legend.title=element_blank(), legend.position = c(0.9,0.85), legend.key=element_blank(), legend.background=element_blank())

png("Figure7A_WCAallresults.png",width=5, height=6,res = 400,units = "in")
plot(p2)
dev.off()

### Things to plot

### Duplicate samples
num_non_miss num_match indiv_1 indiv_2 collection_1 collection_2 sample_type_1 repunit_1 sample_type_2 repunit_2
So, number of unique?
124884 pairs in the data matching

#mixing proportion not super relevant here

#indiv posteriors ## msot important
PofZ  posterior means of likeligood
so where pofZ isn't zero -> highest value
mixture_collection indiv repunit collection PofZ log_likelihood z_score n_non_miss_loci n_miss_loci
Yaounde_Nkolndongo_Market Y200_1 Makokou Makokou 0 -1805.38855744415 37.6113093913832 14088 398

# also for BS

# z scores
If you see a z-score to the maximum-a-posteriori population for an individual in your mixture sample that is considerably less than z_scores you saw in the reference, then you might infer that the individual doesn’t actually fit any of the populations in the reference well.
#but I don't have reference Z scores ? to discuss with Jordi

##Self assigment
#all inferred with likligood for each
#something about most likely ?  or populations with highest self assingment?
indiv collection repunit inferred_repunit repu_scaled_like
12B_1 Mt.Kouffe Mt.Kouffe Abong_Mbang 8.94298459180414e-127

###BONE
#highest rpob for solpath an dWTA
Abong_Mbang Adele Akonolinga Amessodjikope Anguma Asejire_Market Assoukoko Banyo Bayib-Assibong Bayomen Bioko_Batete Bipindi Bisobinam Bongoro Campo Come Dassioko Dimbokro Djoum Ebenguan Ekombitie Emangos Eseka Esse Foumbot Gnanhouizounme Gnkoltang Grand_Lahou Ituri_Forest_RAFALE_Site Kakum_NP Kamina Lolodorf Maan Makokou Manengole Manyemen_Nguti Medouneu Miki Misergue Mt.Kouffe Nditam Oyem Sangmelima SE_Guinea Taguete Tegon_Lama_Forest_Reserve Toumodi Tshuapa_near_Quatorze Uma Yabassi DifferenceProportion
#
outsiders= "DifferenceProportion" is the MSE between probability of the origin estimated with the solution path method and expected probability of the origin. = (top 10 lowest MSE between the estimated and the expected probability of the origin)

#also provides mixing proption

###for now= plot on barplot like ngsadmix
very easy for BONE / then pull out pofZ for rubias
## thne look at MSE