Npa.sim <- function (x, n.sim = 999, N=5) 
{
  require(adegenet)
  require(poppr)
  #if (class(x) != "genelight") #(!is.genind(x)) 
   # stop("x is not a valid genlight object")
  n.pop <- seppop(x)
  f0<-function(z,N){
    selsamp<-sample(z@pop,N) #rownames(z@tab),N)
    z <- z[z@pop %in% selsamp]} #row.names(z@tab)
    #z <- df2genind(as.data.frame(z), sep="/")}
  f1<-function(n.pop,N){
    gd<-n.pop
    for(i in 1:length(n.pop)){
      gd[[i]]<-f0(gd[[i]],N)}
      zgd<-repool(gd)
      #print(zgd)
      x.pA<-private_alleles(zgd)
      x.pA.pop<-c(1:length(x.pA[,1]))
      for(i in 1:length(x.pA[,1])){x.pA.pop[i]<-length(x.pA[i,which(x.pA[i,]>=1)])}
    return(x.pA.pop)}
  f2<-function(x,N,n.sim){
    sim<-t(replicate(n.sim, f1(n.pop,N=N))) 
    return(sim)}
  simf2<-as.data.frame(f2(x,N,n.sim))
  colnames(simf2)<-names(n.pop)
  pop.npa.sim<-stack(simf2)
  colnames(pop.npa.sim)<-c("npa","pop"); 
  return(pop.npa.sim)
}
