#!/bin/sh
#!/usr/bin/env Rscript

#module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R
#Rscript $bin/plot_seed_locator.R $loc_dir/ seed_prediction_results.txt



library(ggplot2)
library(rgdal)
library(plyr)
library(dplyr)
library(raster)
library(ggnewscale)
library(RColorBrewer)
library(ggpubr)

loc_dir <- args[1]
pred_file <- args[2]

print(c(loc_dir, pred_file))

#for non-iteration
#loc_dir <- "/work/idumville/pangolins/RAD/results/locator/seeding.concat.allWCA/"
#pred_file <- "seed_prediction_results.txt"



setwd(loc_dir)

predicted<-read.delim(paste("./", pred_file,sep=""),header=T, fill=T, stringsAsFactors=F)
predicted$opprank <- 1 - as.numeric(predicted$rank)
head(predicted)
predicted00 <- predicted[which(predicted$filter == "00"),]
predicted90 <- predicted[which(predicted$filter == "90"),]
predictedknown <- predicted[which(predicted$coord == "collection"),]



map <- readRDS("/work/idumville/pangolins/RAD/data/geodata/smallmap.rdata")

lon <- c(8.000, 14.2194)
lat <- c(-0.271891, 7.25)

samples <- sample(predicted90$sample, size=1)
samples <- unique(predicted90$sample)

pdf("./locator_individual_plots.pdf")#,onefile = TRUE)
for (indv in samples) {
    predictedtemp <- predicted[which(predicted$sample == indv),]
    predicted00 <- predictedtemp[which(predictedtemp$filter == "00"),] %>%                                      # Top N highest values by group
  arrange(desc(rank)) %>% 
  group_by(sample) %>%
  slice(1:100) 
  predicted90 <- predictedtemp[which(predictedtemp$filter == "90"),] %>%                                      # Top N highest values by group
  arrange(desc(rank)) %>% 
  group_by(sample) %>%
  slice(1:100) 
    predictedknown <- predictedtemp[which(predictedtemp$coord == "collection"),]
p00 <- map + new_scale_fill() + #scale_colour_brewer(palette = "Dark2") +
#scale_fill_brewer(palette = "Greens") + 
geom_point(data=predicted00, aes(x=as.numeric(y), y=as.numeric(x), fill=as.numeric(opprank), colour=as.numeric(opprank)), size=0.6, pch=21) + 
geom_point(data=predictedknown, aes(x=as.numeric(y), y=as.numeric(x)), size=1.2, pch=21) + 
theme(legend.position = "bottom") + theme_bw() + xlim(lon) + ylim(lat)

p90 <- map + new_scale_fill() + #+ scale_colour_brewer(palette = "Greens") +
#scale_fill_brewer(palette = "Dark2") + 
geom_point(data=predicted90, mapping=aes(x=as.numeric(y), y=as.numeric(x), fill=as.numeric(opprank), colour=as.numeric(opprank)), size=0.6, pch=21) + 
geom_point(data=predictedknown, aes(x=as.numeric(y), y=as.numeric(x)), size=1.2, pch=21) + 
theme(legend.position = "bottom") + theme_bw() + xlim(lon) + ylim(lat)

   print(ggarrange(p00, p90, 
          labels = c(paste("No Filtering", indv) , paste("90% loci filtering", indv)), common.legend = TRUE, legend="none",
          ncol = 2, nrow = 1))
}
dev.off()


#making png per indivudal
predicted90 <- predicted[which(predicted$filter == "90"),]
samples <- unique(predicted90$sample)


for (indv in samples) {
png(paste0("./locator_individual_plots_", indv, ".png"))
    predictedtemp <- predicted[which(predicted$sample == indv),]
    predicted00 <- predictedtemp[which(predictedtemp$filter == "00"),] %>%                                      # Top N highest values by group
  arrange(desc(rank)) %>% 
  group_by(sample) %>%
  slice(1:100) 
  predicted90 <- predictedtemp[which(predictedtemp$filter == "90"),] %>%                                      # Top N highest values by group
  arrange(desc(rank)) %>% 
  group_by(sample) %>%
  slice(1:100) 
    predictedknown <- predictedtemp[which(predictedtemp$coord == "collection"),]
p00 <- map + new_scale_fill() + #scale_colour_brewer(palette = "Dark2") +
#scale_fill_brewer(palette = "Greens") + 
geom_point(data=predicted00, aes(x=as.numeric(y), y=as.numeric(x), fill=as.numeric(opprank), colour=as.numeric(opprank)), size=0.6, pch=21) + 
geom_point(data=predictedknown, aes(x=as.numeric(y), y=as.numeric(x)), size=1.2, pch=21) + 
theme(legend.position = "bottom", plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + theme_bw() + xlim(lon) + ylim(lat)

p90 <- map + new_scale_fill() + #+ scale_colour_brewer(palette = "Greens") +
#scale_fill_brewer(palette = "Dark2") + 
geom_point(data=predicted90, mapping=aes(x=as.numeric(y), y=as.numeric(x), fill=as.numeric(opprank), colour=as.numeric(opprank)), size=0.6, pch=21) + 
geom_point(data=predictedknown, aes(x=as.numeric(y), y=as.numeric(x)), size=1.2, pch=21) + 
theme(legend.position = "bottom", plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + theme_bw() + xlim(lon) + ylim(lat)
   print(ggarrange(p00, p90, 
          labels = c(paste("No Filtering", indv) , paste("90% loci filtering", indv)), common.legend = TRUE, legend="none", font.label = list(size = 11), ncol = 2, nrow = 1) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")))
dev.off()
}



#only 90%
pdf("./locator_individual_plots90.pdf")#,onefile = TRUE)
for (indv in sample(samples, size=10)) { #10 random samples
    predictedtemp <- predicted[which(predicted$sample == indv),]
    predicted00 <- predictedtemp[which(predictedtemp$filter == "00"),]
    predicted90 <- predictedtemp[which(predictedtemp$filter == "90"),]
    predictedknown <- predictedtemp[which(predictedtemp$coord == "collection"),]

p90 <- map + new_scale_fill() + #new_scale_colour() + scale_colour_brewer(palette = "Greens") +
#scale_fill_brewer(palette = "Dark2") + 
geom_point(data=predicted90, mapping=aes(x=as.numeric(y), y=as.numeric(x), fill=as.numeric(rank),  colour=as.numeric(rank)), size=0.6, pch=21) + 
#geom_point(data=predictedknown, mapping=aes(x=y, y=x, fill="black"), size=0.8, colour="black",pch=21) + 
theme(legend.position = "bottom") + theme_bw() + xlim(lon) + ylim(lat)

}
dev.off()


