#!/bin/sh
#!/usr/bin/env Rscript

# This is potentially called by loci.sh but no longer available

args = commandArgs(trailingOnly=TRUE)

indv=args[1]
cov_dir=args[2]

setwd(cov_dir)
print(cov_dir)

tab<-read.table(paste0(indv,".f1.bamhits.scaffoldsep"),header=F)

tab[,1]<-as.numeric(tab[,1])
tab[,5]<-as.numeric(tab[,5])
head(tab)

#glengths<-read.table(file="../genomelengths.tsv",header=F)
#tab <- merge(glengths, tab, all=TRUE,by.x="V3",by.y="V4")

if(sum(tab[,1])<4){print("not enough data, quiting"); q(save='no')}

x <- as.data.frame(table(tab[,4]))
#number of scaffolds
m <- length(x$Var1)
tab2 <- tab[order(tab[,4]),]
tab2[,1] <- as.numeric(tab2[,1])
x <- as.data.frame(table(tab2[,4]))


#getting first scaffold with less than 20 loci
x[] <- lapply(x, as.numeric)
p <- min(which(x$Freq < 20))
x2 <-x[which(x$Var1 > p),]
m2 <- length(x2$Freq)
x3 <-x[which(x$Var1 < p),]


#also want another as percentages
x4 <- as.data.frame(table(tab[,4]))
x4[] <- lapply(x4, as.numeric)
x4$prop <- x4$Freq/length(tab[,4])

x5 <-x4[which(x4$Var1 < p),]
x6 <-x4[which(x4$Var1 > p),]

#including length of scaffolds
glengths<-read.table(file="../genomelengths.tsv",header=F)
glengths$V4 <- glengths$V4/1000

tabble <- merge(glengths, x, all.y=TRUE,by.x="V3",by.y="Var1")

tabble$locithou <- tabble$Freq/tabble$V4

tabble$V3<-as.numeric(tabble$V3)
tabble$locithou<-as.numeric(tabble$locithou)

print("scaffold number, with length in kbp, GC%, number of loci, number of loci per kbp")
head(tabble)
ma <- which(tabble$locithou == max(tabble$locithou))
tabble2 <- tabble[-ma,]

vma <- max(tabble$locithou)
#making file
pdf(paste0(indv,".f1.loci.pdf"))
barplot(height=x$Freq, names=x$Var1, space = 0.05,  cex.names=0.5, las=2, xlab="Scaffold", ylab="Number of loci")

#making second table where scaffolds from first one with less than 20 loci are excluded

barplot(height=x2$Freq, names=x2$Var1, space = 0.05, cex.names=0.5, las=2, xlab=paste0("Scaffolds above ", p), ylab="Number of loci",
        main=paste0(indv," scaffolds above number ", p))

barplot(height=x3$Freq, names=x3$Var1, space = 0.05, cex.names=0.5, las=2, xlab=paste0("Scaffolds below ", p), ylab="Number of loci",
        main=paste0(indv," scaffolds below number ", p))

#m <- length(proportions$Freq)
barplot(height=x5$prop, names=x5$Var1, space = 0.05,  cex.names=0.5, las=2, xlab="Scaffold", ylab="Proportion of loci", main = "larger scaffolds")
barplot(height=x6$prop, names=x6$Var1, space = 0.05,  cex.names=0.5, las=2, xlab="Scaffold", ylab="Proportion of loci", main = "smaller scaffolds")
barplot(height=x4$prop, names=x4$Var1, space = 0.05,  cex.names=0.5, las=2, xlab="Scaffold", ylab="Proportion of loci", main = "all scaffolds")

barplot(height=tabble2$locithou, names=tabble2$Var1, space = 0.05,cex.names=0.5, las=2, xlab="Scaffold",  ylab="Number of loci per kbp", main=paste0(indv, "Number of loci per contig wo ", ma))

barplot(height=tabble$locithou, names=tabble$Var1, space = 0.05,cex.names=0.5, las=2, xlab="Scaffold", ylab="Number of loci per kbp")

dev.off()

#also want to print the indv and the p to file
string<-paste0(indv,"_",p,"\n",ma,vma)
write(string, file = "../individuals_scaffold_dropped.txt", append=TRUE)


