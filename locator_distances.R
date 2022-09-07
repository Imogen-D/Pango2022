library(geodist)
library(data.table)

#infile="/work/idumville/pangolins/RAD/results/locatorAL/B51_1_ALminmac3_lociindv/correctpredlocs"

args = commandArgs(trailingOnly=TRUE)
infile=args[1]
indv=args[2]
data=args[3]

setwd("/work/idumville/pangolins/RAD/results/locator")

print("loading data")
if(grepl("predlocs.txt",infile)){
  pd <- fread(infile,data.table=F)
  names(pd) <- c('xpred','ypred','sampleID')
  files <- infile
} else {
  files <- list.files(infile,full.names = T)
  files <- grep("predlocs",files,value=T)
  pd <- fread(files[1],data.table=F)[0,1:3]
  for(f in files){
    a <- fread(f,data.table = F,header=T)[,1:3]
    pd <- rbind(pd,a)
  }
  names(pd) <- c('xpred','ypred','sampleID')
}

alldist <- geodist(pd, measure="geodesic")
av <- mean(alldist)
q <- quantile(alldist, probs = c(.25, .5, .75))
x <- max(alldist)
n <- min(alldist)

cat(indv, av, q, x, n, "\n", sep="\t",file=paste0(data, "locatorsummary.txt"),append=TRUE)