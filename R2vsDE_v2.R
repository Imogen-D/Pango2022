R2vsDE=function(mylist,fmla,DEvarname,lagdist,addplot,col='red',lty=1,minX,maxY){
  
  prep1 = function(x){ x=data.matrix(x)
  x=x[upper.tri(x, diag = FALSE)]}
  
  data=lapply(mylist, prep1) 
  data=as.data.frame(data)
  if (minX<min(data[,colnames(data)==DEvarname],na.rm = TRUE)){
    print(paste0("Minimum distance is ",min(DE)))
    x <- readline("Enter new minimum distance for plot : " )
    minX=as.numeric(unlist(strsplit(x, ",")))
  }
  lag=seq(minX,max(data[,colnames(data)==DEvarname],na.rm = TRUE),lagdist)
  rr=as.data.frame(matrix(0,ncol=2,nrow=length(lag)))
  for (i in c(1:length(lag))){
    datatemp=data
    datatemp2=datatemp[datatemp[,colnames(datatemp)==DEvarname]<lag[i],]
    ca1=summary(lm(as.formula(fmla),datatemp2))
    rr[i,1]=lag[i]
    rr[i,2]=ca1$r.squared
  }
  colnames(rr)=c('Distance','R2')
  if (addplot==FALSE){
    plot(rr,pch=16,type='l',xlab=substitute(paste("Maximum Euclidean distance (",italic('S'),") in Km.",sep="")),
         ylab=expression(italic(R)^2),
         lty=lty,col=col,xlim=c(minX,max(data[,colnames(data)==DEvarname],na.rm = TRUE)),ylim=c(0,maxY))
  }else{
    points(rr,pch=16,type='l',lty=lty,col=col)
  }
  return(rr)
}