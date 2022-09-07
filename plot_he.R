module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3
cd $res/${loc}_indv_he



#for making this into an actual script...

#!/bin/sh
#!/usr/bin/env Rscript
## module load system/R-3.6.1 ; R
# module load system/R-3.5.2 ; module load system/R-3.5.2 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R
## module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R
library(tidyverse)
library(hrbrthemes)
library(viridis)

prefix="rad.ptri"
popmap="/work/idumville/pangolins/RAD/results/autosomes_indv_he/all_indv_he_autosomes.txt"
loc="/work/idumville/pangolins/RAD/results/autosomes_indv_he"
prefix_sample_file="autosomes"
#setwd(loc)


pdf(paste("./",prefix_sample_file,"_he_plot.pdf",sep=""), 50, 10)
#ID/popn/nloci/usedfwread/meancov/qmapreads/lat/long/he/ho
samples<-read.table(popmap,header=F,sep="\t");
samples %>%
  ggplot( aes(x=V2, y=V10, fill=V2)) +
    geom_boxplot() + #violin() +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    #theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none",
      plot.title = element_text(size=11)) +
    ggtitle("Plot of He per individual by population") +
    xlab("Populations") +
    ylab("He per 1000bp")
dev.off()