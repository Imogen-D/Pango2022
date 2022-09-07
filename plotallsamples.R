#!/usr/bin/env Rscript

#bash stuff
awk '{print $25, $6, $7, $39, $32}' /work/jsalmona/pangolins/RAD/data/RAD_pango_recap.txt | sed 's/ /\t/g' > /work/idumville/pangolins/RAD/results/allmarkets.txt

#module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R

library(ggplot2)
library(rgdal)
library(plyr)
library(dplyr)
library(raster)
library(ggnewscale)
library(RColorBrewer)
library(ggrepel)

map <- readRDS("/work/idumville/pangolins/RAD/data/geodata/bigmap.rdata")

setwd("/work/idumville/pangolins/RAD/results/")

#this is super ugly but does the job
meta <- read.table("/work/jsalmona/pangolins/RAD/data/RAD_pango_metadata.txt", header=T)
marks <- meta[,c("Longitude", "Latitude", "Locality", "Market_Type", "lineage")]
ag <- aggregate(.~Locality,  meta[,c("Locality", "new_name")], FUN=head, 1)
count <- meta %>% count(Locality)
ag$count <- count$n
ag <- distinct(merge(ag, marks, by="Locality", all.y=FALSE))

ag$Latiude <- as.numeric(ag$Latitude)
ag$Longitude <- as.numeric(ag$Longitude)

ag$Longitude[which(ag$Market_Type == "Seizure")] <- 3

ag <- ag[-which(ag$Locality == "Banyo"),]

ag$Market_Type.y <- 1
ag$Market_Type.y[which(ag$Market_Type == "Rural_Market")] <- 5

ag$count <- ag$count/1000 #ag$count/1000


ag$Colour <- "DG"
ag$Colour[which(ag$lineage %in% c("WAf", "Gha"))] <- "West"
ag$Colour[which(ag$lineage == c("WCA"))] <- "WCASouth"
ag$Colour[which(ag$lineage == c("Gab"))] <- "Gabon"
ag$Colour[which(ag$lineage == c("CA"))] <- "CA"
ag$Colour[which(ag$Locality %in% c("Bioko_Batete", "Douala_Dakat_Market",  "Douala_Central_Market", "Nditam", "Yabassi",
                                                                   "Bayib-Assibong", "Manyemen_Nguti" , "Foumbot", "Bayomen" ))] <- "WCANorth"
ag$Colour[which(ag$Locality == "Yaounde_Nkolndongo_Market")] <- "WCASouth"

ag <- ag[!duplicated(ag$Locality),]
                                                        
ag$count[which(ag$Locality == "Yaounde_Nkolndongo_Market")] <- ag$count[which(ag$Locality == "Douala_Dakat_Market")]
                                                                   
numsamp <- map + new_scale_fill() +  scale_fill_brewer(palette = "Dark2") + 
scale_color_brewer(palette = "Dark2") + geom_point(ag, mapping=aes(y=Latitude, x=Longitude, fill=as.character(Colour)), shape=21, size = 1.2,  stroke=0.5 ) +
 geom_label_repel(ag, mapping=aes(y=Latitude, x=Longitude, label=new_name, colour = Market_Type), size=2, segment.color="grey20", segment.size=0.5,   label.padding = unit(0.1, "lines"), fontface = "bold", force=0.1, point.size = NA, max.overlaps = Inf)  +
  theme_bw() +  theme(legend.title = element_blank(),
        legend.position=c(x=0.05, y=0.2),
        legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'), #change legend key width
        legend.text = element_text(size=5),
        legend.background=element_blank()) + 
  scale_colour_manual(values = c(Rural_Market = "black", Urban_Market="red", Seizure="orange")) +
  scale_fill_manual(values = c(West =  "yellow2", #"#FFFFB3", 
                                 DG = "#BEBADA",
                                 WCANorth = "#80B1D3", WCASouth =  "#FDB462", Gabon ="#FB8072", CA = "#FB8072")) 
                                 
png("Supp_FigS1_numbersamples.png",width=10, height=4,res = 600,units = "in") 
plot(numsamp)
dev.off()
