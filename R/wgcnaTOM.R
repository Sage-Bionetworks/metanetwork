wgcnaTOM <- function(data,path=NULL,pval=1,outputpath,RsquaredCut=.80,defaultNaPower=6){
  library(WGCNA)
  res <- WGCNA::pickSoftThreshold(data,RsquaredCut=RsquaredCut)
  if(is.na(res$powerEstimate)){
    res$powerEstimate<-defaultNaPower
  }
  network <- abs(cor(data))^res$powerEstimate
  #save(network,file='result_wgcnaSoftThreshold.rda')
  print(res$powerEstimate)
  network <- TOMsimilarity(network)
  #save(network,file=paste0(outputpath,'result_wgcnaTOM.rda'))
  write.csv(network,file=paste0(outputpath,'wgcnaTopologicalOverlapMatrixNetwork.csv'),quote=F)
}
