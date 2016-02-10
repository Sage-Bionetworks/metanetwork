wgcnaTOM <- function(data,path=NULL,pval=1,outputpath,RsquaredCut=.80,defaultNaPower=6){
  library(WGCNA)
  res <- WGCNA::pickSoftThreshold(data,RsquaredCut=RsquaredCut)
  if(is.na(res$powerEstimate)){
    res$powerEstimate<-defaultNaPower
  }
  network <- abs(cor(data))^res$powerEstimate
  #save(network,file='result_wgcnaSoftThreshold.rda')
  print(res$powerEstimate)
  cn <- colnames(network)
  network <- TOMsimilarity(network)
  colnames(network) <- cn
  rownames(network) <- cn
  #save(network,file=paste0(outputpath,'result_wgcnaTOM.rda'))
  network <- network*upper.tri(network)
  write.csv(network,file=paste0(outputpath,'wgcnaTopologicalOverlapMatrixNetwork.csv'),quote=F)
}
