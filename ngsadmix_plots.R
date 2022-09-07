module load system/R-3.5.2 ; module load system/R-3.5.2 libraries/gdal-2.3.0 libraries/proj-4.9.3 ; R

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

#NOW FOR MITO (change dir)

ngsadmix_dir="/work/idumville/pangolins/RAD/results/ngsadmix.mito"
sample_file="pango.rad.all_mito"
order_file="/work/idumville/pangolins/RAD/data/pango_sample_list.txt"
prefix_sample_file="rad.ptri.pango.rad.all_mito"
prefix="rad.ptri"
popmap="/work/idumville/pangolins/RAD/results/ngsadmix.mito/sampledatamito.txt"
loc="/work/idumville/pangolins/RAD/results/ngsadmix.mito"

#for autosomes
ngsadmix_dir="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes"
sample_file="pango.rad.all_autosomes"
order_file="/work/idumville/pangolins/RAD/data/pango_sample_list.txt"
prefix_sample_file="rad.ptri.pango.rad.all_autosomes"
prefix="rad.ptri"
popmap="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/sampledataautosomes.txt"
loc="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes"


#for autosomes lowreads
ngsadmix_dir="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads"
sample_file="pango.rad.all_autosomes"
order_file="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads/highreadsamples.txt"
prefix_sample_file="lowreads.pango.rad.autosomes"
prefix="rad.ptri"
popmap="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads/sampledataautosomes.txt"
loc="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads"
low_list="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/lowurbanreads.txt"


#for xchrom lowreads
ngsadmix_dir="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom.lowreads"
sample_file="pango.rad.all_xchrom"
order_file="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads/highreadsamples.txt"
prefix_sample_file="lowreads.pango.rad.xchrom"
prefix="rad.ptri"
popmap="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom.lowreads/sampledataxchrom.txt"
loc="/work/idumville/pangolins/RAD/results/ngsadmix.xchrom.lowreads"


library(RColorBrewer)
library("dplyr")
#!/bin/sh
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
popmap<-args[1]
prefix_sample_file<-args[2]
loc<-args[3]
order_file<-args[4]
# prefix_sample_file<-"Xrad.mms.rad.all" ; args[1]<-"../../data/RAD_PG2020.data.info" ; loc="/work/pgaubert/monk_seals/RAD/results/ngsadmix/" ;order_file<-"../../data/RAD_PG2020_sample_list.txt"
setwd(loc)


########## likevalues ##################
pdf(paste("./",prefix_sample_file,"_ngsadmix_likevalues.pdf",sep=""), 10, 10)
##
likevalues<-read.table(paste("./",prefix_sample_file,"_likevalues.txt",sep=""))
runame<-prefix_sample_file
## plot likelihood values
par(mfrow=c(2,2))
plot(likevalues[,2],likevalues[,1],type="p", xlab="K", ylab="Likelihood")
boxplot(likevalues[,1]~likevalues[,2], xlab="K", ylab="Likelihood")
plot(likevalues[,2],likevalues[,3],type="p", xlab="K", ylab="Iterations")
plot(likevalues[,3],likevalues[,1],type="p", xlab="Iterations", ylab="Likelihood")
## deltaK procedure
simsum<-likevalues ; colnames(simsum)<-c("LnPr","K","it","rep")
maxK = max(simsum[,2])
Kmeans = tapply(simsum[,1],simsum[,2],mean)
Kvars = tapply(simsum[,1],simsum[,2],var)
Ksds = Kvars^0.5
deltaK = rep(NA,maxK)
for(x in 2:maxK){deltaK[x] = abs(Kmeans[x+1] - 2*Kmeans[x] + Kmeans[x-1])/Ksds[x]}
par(mfrow=c(2,2))
bp<-boxplot(simsum[,1]~simsum[,2],xlab="K", ylab='Ln')
plot(deltaK,type = 'b',xlab='K',cex.axis=0.75)
#bp<-boxplot(simsum[,1]~simsum[,2],xlab="K", ylab='Ln',log="y")
plot(deltaK,type = 'b',xlab='K',cex.axis=0.75,log="y")
dev.off()


png(paste("./",prefix_sample_file,"_ngsadmix_likevalues_final.png",sep=""), width = 14, height = 8, units = 'in', res = 300)
par(mfrow=c(2,1))
boxplot(likevalues[,1]~likevalues[,2], xlab="K", ylab="Likelihood")
mtext("  A",side=3,line=1,adj=0,font = 2)
plot(deltaK,type = 'b',xlab='K',cex.axis=0.75,log="y")
mtext("  B",side=3,line=1,adj=0,font = 2)
dev.off()

png(paste("./Supp_FigS7_",prefix_sample_file,"_ngsadmix_reads_coverage.png",sep=""), 10, 20, units="in", res=500)
par(mfrow=c(1,2)) #no need for structuring

#pdf(paste("./",prefix_sample_file,"_ngsadmix_reads_assign.pdf",sep=""), 10, 10)
#par(mfrow=c(2,2)) no need for structuring
for(K in 2:20){
i=K
  samples<-read.table(popmap,header=F,sep="\t"); #samples<-samples[-which(samples$V1 == "Y199_1"),] #only for xchrom
  dim2<-(dim(samples)[2])+1
  for(j in 1:length(samples[,1])){ samples[j,dim2]<-paste(samples[j,1],samples[j,4],sep=".")} ; samples[,1]<-samples[,dim(samples)[2]]
  
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
  dim2<-(dim(samples)[2])+1 ; 
  for(j in 1:length(samples[,1])){ samples[j,dim2]<-max(admix[,j])} 
plot(as.numeric(samples$V4),samples$V14,xlab="Number of reads",ylab="NGSadmix assignment",log="x",main=paste("RAD K =",K),ylim=c(0.3,1))
}
#dev.off()


## plot assignment as a function of coverage
#pdf(paste("./",prefix_sample_file,"_ngsadmix_reads_coverage.pdf",sep=""), 10, 10)

for(K in 2:20){
i=K
  samples<-read.table(popmap,header=F,sep="\t"); #samples<-samples[-which(samples$V1 == "Y199_1"),] #only for xchrom
  dim2<-(dim(samples)[2])+1
  for(j in 1:length(samples[,1])){ samples[j,dim2]<-paste(samples[j,1],samples[j,5],sep=".")} ; samples[,1]<-samples[,dim(samples)[2]]
  
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
  dim2<-(dim(samples)[2])+1 ; 
  for(j in 1:length(samples[,1])){ samples[j,dim2]<-max(admix[,j])} 

plot(as.numeric(samples$V5),samples$V14,xlab="Coverage",ylab="NGSadmix assignment",log="x",main=paste("RAD K =",K),ylim=c(0.3,1))
}
dev.off()



#plotting all liklihood values for K=6
dfemp <- data.frame()
library(stringr)
library(tidyr)
for(L in 1:20){
K=6
i=K
  samples<-read.table(popmap,header=F,sep="\t"); #samples<-samples[-which(samples$V1 == "Y199_1"),] #only for xchrom
  dim2<-(dim(samples)[2])+1
  for(j in 1:length(samples[,1])){ samples[j,dim2]<-paste(samples[j,1],samples[j,4],sep=".")} ; samples[,1]<-samples[,dim(samples)[2]]   
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",L,".qopt",sep=""))))
  dim2<-(dim(samples)[2])+1 ; 
  for(j in 1:length(samples[,1])){ samples[j,dim2]<-max(admix[,j])} 
  df <- samples[,c("V13", "V4", "V5", "V14")]
  df$L <- L
  dfemp <- rbind(dfemp,df)
}
dfemp <- separate(dfemp,V13, c("indv", NA), sep="[.]")
data_msd <- dfemp %>%     
  mutate(indv = as.factor(indv))   %>%                 # Get mean & standard deviation by group
  group_by(indv) %>%
  summarise_at(vars(V14),
               list(mean = mean,
                    sd = sd, min=min, max=max)) %>% 
  as.data.frame()
head(data_msd)     

#V5 for coverage, V4 for reads
#pdf(paste("./",prefix_sample_file,"_ngsadmix_reads_assign_all_runs.pdf",sep=""), 18, 10)

png(paste("./Supp_FigS7_",prefix_sample_file,"_ngsadmix_reads_coverage.png",sep=""), 12, 15, units="in", res=500)
par(mfrow=c(2,1)) 
sampletemps <- separate(samples,V13, c("indv", NA), sep="[.]")
data_msd <- merge(data_msd, sampletemps, by="indv")
low_high <- read.delim("/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/Low_High.txt")
data_msd <- merge(data_msd, low_high, by.x="indv", by.y="new_name")
colors= rev(brewer.pal(n = 2, name = "Set2"))
plot(as.numeric(data_msd$V4/1000),data_msd$mean,xlab="Number of Used Reads (x1000)",ylab="NGSadmix assignment",ylim=c(0.2,1.01), log="x", bg = colors[ unclass(data_msd$LowRead) ], pch=21) 
arrows(as.numeric(data_msd$V4/1000),(data_msd$min),as.numeric(data_msd$V4/1000),(data_msd$max), code=3, length=0.02, angle = 90)
abline(v=153, col="red")
legend("topleft", legend = c("High read",  "Low Read"), fill = colors, bty = "n")
mtext("A", 2, adj=1, las=1, line=1, padj=-23, font=2)
plot(as.numeric(data_msd$V5),data_msd$mean,xlab="Mean Coverage",ylab="NGSadmix assignment",ylim=c(0.2,1.01), log="x", bg = colors[ unclass(data_msd$LowRead) ], pch=21) 
arrows(as.numeric(data_msd$V5),(data_msd$min),as.numeric(data_msd$V5),(data_msd$max), code=3, length=0.02, angle = 90)
mtext("B", 2, adj=1, las=1, line=1, padj=-23, font=2)
dev.off()

#plotting liklihood of assingment over binned values
pdf(paste("./",prefix_sample_file,"_ngsadmix_means_reads_assign.pdf",sep=""), 10, 10)
data_msd_binned <- data_msd %>% mutate(x_bins=cut(V4, breaks=20))
ggplot(data=data_msd_binned, aes(x=x_bins, y=mean)) + geom_boxplot() + theme_bw() + xlab("Read Bins") + ylab("Assigment Probability")
data_msd_binned <- data_msd %>% mutate(x_bins=cut(V5, breaks=20))
ggplot(data=data_msd_binned, aes(x=x_bins, y=mean)) + geom_boxplot() + theme_bw() + xlab("Coverage Bins") + ylab("Assigment Probability")
ggplot(data=data_msd, aes(x=V5, y=mean)) + geom_point() + geom_smooth() + scale_x_continuous(trans='log2') + theme_bw() + xlab("Coverage Bins") + ylab("Assigment Probability")
dev.off()


  #want all values



#main=paste("RAD K = 6",L),
#read log="x", and divide by 1000

#gettiing samples first 95% cutoff in reads and coverage
K=6
i=K
  samples<-read.table(popmap,header=F,sep="\t"); #samples<-samples[-which(samples$V1 == "Y199_1"),] #only for xchrom
  dim2<-(dim(samples)[2])+1
  for(j in 1:length(samples[,1])){ samples[j,dim2]<-paste(samples[j,1],samples[j,4],sep=".")} ;# samples[,1]<-samples[,dim(samples)[2]]
  
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
  dim2<-(dim(samples)[2])+1 ; 
  for(j in 1:length(samples[,1])){ samples[j,dim2]<-max(admix[,j])} 
  #sort by reads
  readsort <- samples[order(samples$V4),]
  #readsort[c(1:100), c(1,4,14)]
 #manually looked for first value > 0.95  number reads = 151742 
 lowreads <- as.character(readsort[c(1:25),1])
 temp <- readsort[c(26:548),]
  temp <- as.character(temp[which(temp$V9 %in% c("DG", "Gha", "WAf", "CA")),1])
write.table(c(lowreads, temp), file="./lowreadswoFGGhaWAfCA.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

#Making file for IBD
  temp <- c(lowreads, temp)
  #temp <- samples[which(!samples$V1 %in% temp),]
  Bio <- as.character(samples$V1[which(samples$V2 == "Bioko_Batete")])
  unique(c(temp, Bio))
  Urban <- as.character(samples$V1[which(samples$V10 == "Urban_Market")])
write.table(unique(c(temp, Bio)), file="./lowreadswoFGGhaWAfCAwoBioko.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(unique(c(temp, Bio, Urban)), file="./lowreadswoFGGhaWAfCAwoBiokowoUrban.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)


  #seperately sort by coverage
  #covsort <- samples[order(samples$V5),] #manually looked for first value > 0.95 >1.5
  #covrem <- covsort[c(1:19),]
  #which(covrem$V1 %in% readrem$V1) #complete overlap so just remove read samples

## plot ngsadmix results
pdf(paste("./",prefix_sample_file,"_ngsadmix_nice_barplot.pdf",sep=""), 100,30)
layout(matrix(c(1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19), nrow = 1, ncol = 20, byrow = TRUE)) #layout(matrix(c(1,1,2,3,4,5,6,7,8,9), nrow = 1, ncol = 10, byrow = TRUE))
for(K in 2:20){
i=K
  samples<-read.table(popmap,header=F); #sort(samples[,1])
  #samples<-samples[-which(samples$V1 == "Y199_1"),] #removing y199_1 for xchrom Re-add IN AUTOSOMES
  pg2020<-read.table(order_file,header=F)  
  #pg2020<-pg2020[-which(pg2020$V1 == "Y199_1"),]
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
      poplist <- c("SeizureP", "SeizureB", "Gnkoltang", "Oyem", "Makokou", "GMedouneu", "Bisobinam", "Taguete","Misergue", "Emangos", "Anguma", "Ebenguan", "Bongoro", "Campo", "Maan", "Djoum", "Sangmelima", "Lolodorf", "Bipindi", "Ekombitie", "Eseka", "Akonolinga", "Yaounde_Nkolndongo_Market", "Esse", "Abong_Mbang",  "Douala_Central_Market", "Douala_Dakat_Market", "Yabassi", "Bayomen", "Manengole", "Foumbot", "Bayib-Assibong", "Nditam", "Manyemen_Nguti", "Banyo", "Bioko_Batete", "Mankessim",  "Kumasi")

  
  #poplist <- poplist[which(poplist %in% samples$V2)]
  samples[,2]<-factor(samples[,2],levels=poplist) #; samples<-samples[order(samples[,2]),] # ; rownames(samples) <- NULL
  #samples<-samples[,order(pg2020)]
  #admix<-admix[,order(pg2020)]; #[,1]
  #        admix<-admix[,order(as.numeric(rownames(samples)))]
  admix<-admix[,order(samples[,2])]; 
  samples<-samples[order(samples[,2]),]
          #rownames(samples) <- NULL
  niveaux<-levels(as.factor(samples[,2]))
  pop<-samples[order(match(samples[,2],niveaux)),]
  tempo<-tapply(1:nrow(pop),pop[,2],mean);
  tempo<-tempo[order(tempo)]
  tempomax<-tapply(1:nrow(pop),pop[,2],max);
  tempomax<-tempomax[order(tempomax)]
  


  if(K==2){
    par(mar=c(0.1,5,1.1,0.2),oma=c(0.1,1,0.1,0.1))
    barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", 
            horiz=T,axes=F,names.arg=samples[,1],las=2,cex.names=0.5)
    title(paste("K =",K),line=-2)
    text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=1.3,col="black")
    abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")
  } else if (K==5) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samples$V5, cex.names=1) #names.arg=samples$coverage
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")    
  } else if (K==7) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samples$V1,las=2,cex.names=0.5) #names.arg=samples$name
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")   
  } else if (K==8) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samples$V9,las=2,cex.names=0.5) #names.arg=samples$lineage
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22") 
  } else if (K==13) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samples$V1,las=2,cex.names=0.5) #names.arg=samples$lineage
      title(paste("K =",K),line=-2)
      text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=1.3,col="black")
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22") 
  } else {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              names.arg=rep("",length(samples[,1])))
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")     
  }}
dev.off()

##plotting just rural and high cov

## plot ngsadmix results
pdf(paste("./",prefix_sample_file,"_ngsadmix_rural_barplot.pdf",sep=""), 60, 30) #60,30 for pdf in inch; 6000,3000 pixels png
layout(matrix(c(1,1,2,3,4,5,6,7,8,9,10,11,12,13), nrow = 1, ncol = 12, byrow = TRUE)) #layout(matrix(c(1,1,2,3,4,5,6,7,8,9), nrow = 1, ncol = 10, byrow = TRUE))
for(K in 2:12){
i=K
  samples<-read.table(popmap,header=F); #sort(samples[,1])
  #samples<-samples[-which(samples$V1 == "Y199_1"),] #removing y199_1 for xchrom Re-add IN AUTOSOMES
  pg2020<-read.table(order_file,header=F)  
  #pg2020<-pg2020[-which(pg2020$V1 == "Y199_1"),]
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
      poplist <- c("SeizureP", "SeizureB", "Gnkoltang", "Oyem", "Makokou", "GMedouneu", "Bisobinam", "Taguete","Misergue", "Emangos", "Anguma", "Ebenguan", "Bongoro", "Campo", "Maan", "Djoum", "Sangmelima", "Lolodorf", "Bipindi", "Ekombitie", "Eseka", "Akonolinga", "Yaounde_Nkolndongo_Market", "Esse", "Abong_Mbang",  "Douala_Central_Market", "Douala_Dakat_Market", "Yabassi", "Bayomen", "Manengole", "Foumbot", "Bayib-Assibong", "Nditam", "Manyemen_Nguti", "Banyo", "Bioko_Batete", "Mankessim",  "Kumasi")

  samples[,2]<-factor(samples[,2],levels=poplist) 
  admix<-admix[,order(samples[,2])]; 
  samples<-samples[order(samples[,2]),]
          #rownames(samples) <- NULL
  niveaux<-levels(as.factor(samples[,2]))
  pop<-samples[order(match(samples[,2],niveaux)),]
  
  #HERE REMOVE INDIVIDUALS FROM POP, SAMPLES AND ADMIX FOR RURAL HIGH READS 
  lines <- readLines(low_list)
  samplenum <- which(samples$V1 %in% lines)
  samples <- samples[-c(samplenum),]
  pop <- pop[-c(samplenum),]
  admix <- admix[,-c(samplenum)]

  tempo<-tapply(1:nrow(pop),pop[,2],mean);
  tempo<-tempo[order(tempo)]
  tempomax<-tapply(1:nrow(pop),pop[,2],max);
  tempomax<-tempomax[order(tempomax)]
  

  if(K==2){
    par(mar=c(0.1,5,1.1,0.2),oma=c(0.1,1,0.1,0.1))
    barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", 
            horiz=T,axes=F,names.arg=samples[,1],las=2,cex.names=0.5)
    title(paste("K =",K),line=-2)
    text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=1.3,col="black")
    abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")
  } else if (K==6) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samples$V5, cex.names=1) #names.arg=samples$coverage
      title(paste("K =",K),line=-2)
    text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=1.5,col="black")
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")    
  } else if (K==7) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samples$V1,las=2,cex.names=1) #names.arg=samples$name
      title(paste("K =",K),line=-2)
      text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=1.3,col="black")
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")   
  } else if (K==8) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samples$V9,las=2,cex.names=0.5) #names.arg=samples$lineage
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22") 
  } else {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              names.arg=rep("",length(samples[,1])))
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")     
  }}
dev.off()

##########
  
pdf(paste("./",prefix_sample_file,"_ngsadmix_highcov.pdf",sep=""), 30,30)
layout(matrix(c(1,1,2,3,4,5,6,7,8,9), nrow = 1, ncol = 10, byrow = TRUE))
for(K in 2:10){
i=K
  samples<-read.table(popmap,header=F); sort(samples[,1])
  #ONLY NEED THESE THREE COMMANDS FOR XCHROM WEIRDNESS (MIGHT HAVE TO REMOVE SOME FOR MITO TOO?)
    #samples<-samples[-which(samples$V1 == "Y199_1"),]
    #samples<-droplevels(samples)
    #rownames(samples) <- NULL
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
  samplesb<-samples[which(as.numeric(samples$V5) > 3),] #COVERAGE ABOVE 3
  #FOR TESTING samplesb<-samples[which(as.numeric(samples$V5) > 15),] #coverage above 15
   poplist<-c("Guinea", "Cote_dIvoire", "Ghana", "Togo", "Benin", "Nigeria", "Banyo", "Foumbot", "Nditam", "Manengole", "SO", "Yabassi", "DCM", "DDM", "Bioko", "Eseka", "Bayomen", "Yaounde", "Tongue", "Ekombitie", "Akonolinga", "Abgab", "Abongmbang", "Sangmelima",  "PNCM", "GES", "NorthG", "GFV", "Central_African_Republic", "DR_Congo", "SeizureB", "SeizureP")
  samplesb[,2]<-factor(samplesb[,2],levels=poplist);
  samplesb<-samplesb[order(samplesb[,2]),]
  tokeep <- rownames(samplesb)
  tokeep<-as.numeric(tokeep)
  admix <- admix[,tokeep]
  admix<-admix[,order(samplesb[,2])]
  niveaux<-levels(as.factor(samplesb[,2]))
  pop<-samplesb[order(match(samplesb$V2,niveaux)),] #pop<-samplesb[order(match(samplesb[,2],niveaux)),]
  tempo<-tapply(1:nrow(pop),pop[,2],mean);
  tempo<-tempo[order(tempo)]
  tempomax<-tapply(1:nrow(pop),pop[,2],max);
  tempomax<-tempomax[order(tempomax)]
  if(K==2){
    par(mar=c(0.1,5,1.1,0.2),oma=c(0.1,1,0.1,0.1))
    barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", 
            horiz=T,axes=F,names.arg=samplesb[,1],las=2,cex.names=0.5)
    title(paste("K =",K,"\n","covmedian above 1000"),line=-2)
    text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=1.3,col="black")
    abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")
  } else if (K==5) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samplesb$V5, cex.names=0.5) #names.arg=samples$coverage
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")   
  } else if (K==7) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samplesb$V1,las=2,cex.names=0.5) #names.arg=samples$name
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")   
  } else if (K==8) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samplesb$V9,las=2,cex.names=0.5) #names.arg=samples$lineage
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")     
  } else {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              names.arg=rep("",length(samplesb[,1])))
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")     
  }}
dev.off()



########## Final plot for presentation with Yaounde and Doala
  
ngsadmix_dir="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads"
sample_file="pango.rad.all_autosomes"
order_file="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads/highreadsamples.txt"
prefix_sample_file="lowreads.pango.rad.autosomes"
prefix="rad.ptri"
popmap="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads/sampledataautosomes.txt"
loc="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes.lowreads"
low_list="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/lowurbanreads.txt"


library(RColorBrewer)
library("dplyr")
library(data.table)

setwd(loc)


########## likevalues ##################
likevalues<-read.table(paste("./",prefix_sample_file,"_likevalues.txt",sep=""))
runame<-prefix_sample_file

K=6
i=K
  samples<-read.table(popmap,header=F); sort(samples[,1])
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
  admix <- admix[,c(which(samples$V2 == "Yaounde_Nkolndongo_Market"))]
  samples <- samples[c(which(samples$V2 == "Yaounde_Nkolndongo_Market")),]
    samples<-droplevels(samples)
    rownames(samples) <- NULL
admix <- as.data.frame(admix)[sapply(as.data.frame(admix), function(x) max(x, na.rm = T) > 0.75)]
samples <- samples[c(as.numeric(stringr::str_remove(colnames(admix), "V"))),]

png(paste("./",prefix_sample_file,"_ngsadmix_Yaounde_urban.png",sep=""), 2, 10, units="in", res=600)
layout(matrix(c(1,1), nrow = 1, ncol = 2, byrow = TRUE))
par(mar=c(0.1,5,1.1,0.2),oma=c(0.1,1,0.1,0.1))
    barplot(as.matrix(admix),col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", 
            horiz=T,axes=F,names.arg=samples[,1],las=2,cex.names=0.5)
    #title(paste("K =",K),line=-2)
   # text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=1.3,col="black")
    #abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")  
dev.off()

assign <- stringr::str_remove(colnames(transpose(admix))[max.col(transpose(admix),ties.method="first")], "V")
market <- rep("yao", length(assign))
yao <- data.frame(assign, market, new_name = samples$V1)

samples$V1[which(yao$assign == "3")] #unassinged smaples?
# Y58_1  Y68_1  Y80_1  Y89_1  Y58_2  Y68_2  Y80_2  Y112_1 Y89_2

admix[,which(samples$V1 %in% c("Y199_1" ,"Y360_1", "L20_1" , "L26_1" , "Y66_2" , "Y69_2" , "L27_1" , "LG10_1"))] #BONE no assign

### for douala
K=6
i=K
  samples<-read.table(popmap,header=F); sort(samples[,1])
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
  admix <- admix[,c(which(samples$V2 %in% c("Douala_Central_Market", "Douala_Dakat_Market")))]
  samples <- samples[c(which(samples$V2 %in% c("Douala_Central_Market", "Douala_Dakat_Market"))),]
    samples<-droplevels(samples)
    rownames(samples) <- NULL
admix <- as.data.frame(admix)[sapply(as.data.frame(admix), function(x) max(x, na.rm = T) > 0.75)]
samples <- samples[c(as.numeric(stringr::str_remove(colnames(admix), "V"))),]

png(paste("./",prefix_sample_file,"_ngsadmix_Douala_urban.png",sep=""), 2, 10, units="in", res=600)
layout(matrix(c(1,1), nrow = 1, ncol = 2, byrow = TRUE))
par(mar=c(0.1,5,1.1,0.2),oma=c(0.1,1,0.1,0.1))
    barplot(as.matrix(admix),col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", 
            horiz=T,axes=F,names.arg=samples[,1],las=2,cex.names=0.5)
    #title(paste("K =",K),line=-2)
   # text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=1.3,col="black")
    #abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")  
dev.off()


assign <- stringr::str_remove(colnames(transpose(admix))[max.col(transpose(admix),ties.method="first")], "V")
market <- rep("dou", length(assign))
dou <- data.frame(assign, market, new_name=samples$V1)

write.table(rbind(yao, dou), file="assignmentofurban.txt", row.names = FALSE, quote = FALSE)



################################
# plotting AL Barplot ngsadmix #
################################

likevalues<-read.table(paste("./",prefix_sample_file,"_likevalues.txt",sep=""))
runame<-prefix_sample_file
low_list="/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/ALbarplotsamps.txt" #only rural, removed sangmelima, djoum, foumbot (3 largest)
unique(readLines("/work/idumville/pangolins/RAD/results/ngsadmix.autosomes/ALbarplotlocal.txt")) #only rural, removed sangmelima, djoum, foumbot (3 largest)


#pdf(paste("./",prefix_sample_file,"_ngsadmix_rural_barplot.pdf",sep=""), 60, 30) #60,30 for pdf in inch; 6000,3000 pixels png

png(paste("./Supp_FigS11",prefix_sample_file,"_ngsadmix_AL_barplot.png",sep=""), 10, 10, units="in", res=600) 

layout(matrix(c(1,1,2,3,4,5,6,7,8,9), nrow = 1, ncol = 10, byrow = TRUE))
for(K in 2:10){
i=K
  samples<-read.table(popmap,header=F); 
  pg2020<-read.table(order_file,header=F)  
  #pg2020<-pg2020[-which(pg2020$V1 == "Y199_1"),]
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))

      poplist <- c("Miki", "Uma", "Ituri_Forest_RAFALE_Site", "Tshuapa_near_Quatorze", "Oyem", "Gnkoltang", "Makokou","GMedouneu" , "Taguete" , "Bisobinam",  "Bongoro" , "Ebenguan"  , "Maan", "Misergue",  "Anguma",  "Emangos","Campo",  "Ekombitie", "Abong_Mbang", "Akonolinga" , "Eseka" , "Bipindi" , "Lolodorf", "Manengole", "Esse", "Yabassi" , "Bayomen" , "Manyemen_Nguti", "Nditam" , "Bayib-Assibong", "Bioko_Batete", "Asejire_Market", "Gnanhouizounme" ,"Tegon_Lama_Forest_Reserve", "Mt.Kouffe",  "Come" , "Kamina" , "Amessodjikope",  "Assoukoko", "Adele", "Mankessim", "Kumasi",
 "Kakum_NP" , "Grand_Lahou", "Dimbokro", "Dassioko" , "Toumodi", "SE_Guinea" )
          

  samples[,2]<-factor(samples[,2],levels=poplist) 
  admix<-admix[,order(samples[,2])]; 
  samples<-samples[order(samples[,2]),]
  niveaux<-levels(as.factor(samples[,2]))
  pop<-samples[order(match(samples[,2],niveaux)),]
  
  #HERE KEEP INDIVIDUALS FROM RURAL ONLY
  lines <- readLines(low_list)
  samplenum <- which(samples$V1 %in% lines)
  samples <- samples[c(samplenum),]
  pop <- pop[c(samplenum),]
  admix <- admix[,c(samplenum)]

  tempo<-tapply(1:nrow(pop),pop[,2],mean);
  tempo<-tempo[order(tempo)]
  tempomax<-tapply(1:nrow(pop),pop[,2],max);
  tempomax<-tempomax[order(tempomax)]
  

  if(K==2){
    par(mar=c(0.1,5,1.1,0.2),oma=c(0.1,1,0.1,0.1))
    barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", 
            horiz=T,axes=F,names.arg=samples[,1],las=2,cex.names=0.5)
    title(paste("K =",K),line=-1, cex = 1)
    text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=0.6,col="black")
    abline(h=c(0,tempomax),lty=1,lwd=1,col="grey22")
  } else if (K==5) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samples$V9, cex.names=0.5) #names.arg=samples$coverage
      title(paste("K =",K),line=-1, cex = 1)
    text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=0.6,col="black")
      abline(h=c(0,tempomax),lty=1,lwd=1,col="grey22")    
  } else if (K==8) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              las=2, names.arg=samples$V1,las=2,cex.names=0.5) #names.arg=samples$name
      title(paste("K =",K),line=-1, cex = 1)
      text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=0.6,col="black")
      abline(h=c(0,tempomax),lty=1,lwd=1,col="grey22")   
  } else {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              names.arg=rep("",length(samples[,1])))
      title(paste("K =",K),line=-1, cex = 1)
      abline(h=c(0,tempomax),lty=1,lwd=1,col="grey22")     
  }}
dev.off()



#SCRAP onwards, DONT USE PG2020 ANYMORE
#pg2020<-read.table(order_file,header=F)  
  #pg2020<-pg2020[-which(pg2020$V1 == "Y199_1"),]
   #admix<-admix[,order(pg2020)]; #[,1]
  #admix<-admix[,order(samples[,2])]; 
##
pdf(paste("./",prefix_sample_file,"_ngsadmix_nice_barplot_no-low-cov.pdf",sep=""), 10, 10)
layout(matrix(c(1,1,2,3,4,5,6,7), nrow = 1, ncol = 8, byrow = TRUE))
for(K in 2:8){
i=K
  samples<-read.table(popmap,header=T); sort(samples[,3])
  pg2020<-read.table(order_file,header=F)  
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
  poplist<-c("Abgab", "Abongmbang", "Akonolinga", "Banyo", "Bayomen", "Benin", "Bioko", "Central_African_Republic", "Cote_dIvoire", "DCM", "DDM", "DR_Congo", "Ekombitie", "Eseka", "Foumbot", "GES", "GFV", "Ghana", "Guinea", "Manengole", "Nditam", "Nigeria", "NorthG", "PNCM", "Sangmelima", "SeizureB", "SeizureP", "SO", "Togo", "Tongue", "Yabassi", "Yaounde")
  samples[,2]<-factor(samples[,2],levels=poplist)
  admix<-admix[,order(pg2020[,1])]; 
  admix<-admix[,order(samples[,2])]; 
  samples<-samples[order(samples[,2]),]
    admix<-admix[,samples$coverage>2]
    samples<-samples[samples$coverage>2,]
    admix<-admix[,samples$ID!="MS482"]
    samples<-samples[samples$ID!="MS482",]
  niveaux<-levels(as.factor(samples[,4]))
  pop<-samples[order(match(samples[,4],niveaux)),]
  tempo<-tapply(1:nrow(pop),pop[,4],mean);
  tempo<-tempo[order(tempo)]
  tempomax<-tapply(1:nrow(pop),pop[,4],max);
  tempomax<-tempomax[order(tempomax)]
  if(K==2){
    par(mar=c(0.1,5,1.1,0.2),oma=c(0.1,1,0.1,0.1))
    barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", 
            horiz=T,axes=F ,names.arg=samples[,1],las=2,cex.names=1)
    title(paste("K =",K),line=-2)
    text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=1.3,col="black")
    abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")
  } else if (K==5) {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              names.arg=samples$coverage,las=2,cex.names=1)
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")     
  } else {
    par(mar=c(0.1,0.1,1.1,0.1))
      barplot(admix,col=rainbow(K),space=0,border=NA,ylab="",main="", horiz=T,axes=F ,
              names.arg=rep("",length(samples[,1])))
      title(paste("K =",K),line=-2)
      abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")     
  }}
dev.off()

q(save="no")


pdf(paste("./",prefix_sample_file,"_ngsadmix_barplot.pdf",sep=""), 10, 10)
par(mfrow=c(4,1))
for(i in 2:max(likevalues[,2])){
  samples<-read.table(popmap,header=T); sort(samples[,3])
  pg2020<-read.table(order_file,header=F)  
  bestlike<-which(likevalues[likevalues[,2]==i,1]==max(likevalues[likevalues[,2]==i,1])) #select the best likelihood of the N runs of     
  admix<-t(as.matrix(read.table(paste("./",prefix_sample_file,"_k",i,"/",prefix_sample_file,"_k",i,"_s",bestlike,".qopt",sep=""))))
  poplist<-c("Abgab", "Abongmbang", "Akonolinga", "Banyo", "Bayomen", "Benin", "Bioko", "Central_African_Republic", "Cote_dIvoire", "DCM", "DDM", "DR_Congo", "Ekombitie", "Eseka", "Foumbot", "GES", "GFV", "Ghana", "Guinea", "Manengole", "Nditam", "Nigeria", "NorthG", "PNCM", "Sangmelima", "SeizureB", "SeizureP", "SO", "Togo", "Tongue", "Yabassi", "Yaounde")
  samples[,2]<-factor(samples[,2],levels=poplist)
  admix<-admix[,order(pg2020[,1])]; 
  admix<-admix[,order(samples[,2])]; 
  samples<-samples[order(samples[,2]),]
  niveaux<-levels(as.factor(samples[,2]))
  pop<-samples[order(match(samples[,2],niveaux)),]
  tempo<-tapply(1:nrow(pop),pop[,2],mean);
  tempo<-tempo[order(tempo)]
  tempomax<-tapply(1:nrow(pop),pop[,2],max);
  tempomax<-tempomax[order(tempomax)]
     par(mar=c(0.1,5,1.1,0.2),oma=c(0.1,1,0.1,0.1))
    barplot(admix,col=brewer.pal(n = K, name = "Set3"),space=0,border=NA,ylab="",main="", 
            horiz=T,axes=F,names.arg=samples[,1],las=2,cex.names=0.5)
    title(paste("K =",K),line=-2)
    text(0.5,tempo-0.5,labels=names(tempo),xpd=T,srt=0,cex=1.3,col="black")
    abline(h=c(0,tempomax),lty=1,lwd=1.5,col="grey22")}
dev.off()

