wgcnaSoftThreshold <- function(data,path=NULL,outputpath,RsquaredCut=.80,defaultNaPower=6){
  library(WGCNA)
  res <- WGCNA::pickSoftThreshold(data,RsquaredCut=RsquaredCut)
  if(is.na(res$powerEstimate)){
    res$powerEstimate<-defaultNaPower
  }
  
  network <- abs(cor(data))^res$powerEstimate

  save(network,file=paste0(outputpath,'result_wgcnaST.rda'))
}
