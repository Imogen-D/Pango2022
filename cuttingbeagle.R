
#!/bin/sh
#!/usr/bin/env Rscript


#module load system/R-3.4.3 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R

#library(tidyverse)
library(R.utils)


args = commandArgs(trailingOnly=TRUE)

loc=args[1] #xchrom
fileprefix=args[2] #ruralhighread // name of output beagle
low_list=args[3] #"/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/lowurbanreads.txt"
beagle=args[4] #"/work/idumville/pangolins/RAD/results/ngsadmix.xchrom.lowreads/unfilteredallxchrom.beagle.gz"
popmap=args[5] #"/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/sampledataautosomes.txt"
ngsadmix_dir=args[6] #"/work/idumville/pangolins/RAD/results/ngsadmix.xchrom.lowreads/"

message("arguemnts loaded")
message(c(loc, fileprefix, low_list, beagle, popmap, ngsadmix_dir))
message(paste0(c(ngsadmix_dir, loc, fileprefix, ".beagle.gz")))

setwd(ngsadmix_dir)

samples<-read.table(popmap,header=F)
beagletab <- read.table(gzfile(beagle))  
lines <- readLines(low_list)


samples<-read.table(popmap,header=F)
beagletab <- read.table(gzfile(beagle))  
lines <- readLines(low_list)

samples<-read.table(popmap,header=F)
lines <- readLines(low_list)
samplenum <- which(samples$V1 %in% lines)
x <- samples$V1[-c(samplenum)]
write.table(x, paste0(loc, fileprefix, "samples", ".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)

q(save="no")


if (loc=="xchrom") {
samples <- samples[which(samples$V1 != "Y199_1"),]
}

samplenum <- which(samples$V1 %in% lines)
samplenum <- samplenum - 1

samplenum <- paste("Ind", samplenum, sep="")
message(samplenum)
samplenum <- which(unlist(beagletab[1,], use.names=FALSE) %in% noquote(samplenum))

beagletab2 <- beagletab[,-c(samplenum)]

head(beagletab2)
filename=paste0(loc, fileprefix, "beagle.gz")
message(filename)

write_delim(beagletab2, filename, col_names=FALSE, delim="\t")


samples<-read.table(popmap,header=F)
beagletab <- read.table(gzfile(beagle))  
lines <- readLines(low_list)

samples<-read.table(popmap,header=F)
lines <- readLines(low_list)
samplenum <- which(samples$V1 %in% lines)
x <- samples$V1[-c(samplenum)]
write.table(x, paste0(loc, fileprefix, "samples", ".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)

q(save="no")

