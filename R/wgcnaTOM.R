wgcnaTOM <- function(data,path=NULL,outputpath,RsquaredCut=.80){
  library(WGCNA)
  res <- WGCNA::pickSoftThreshold(data,RsquaredCut=RsquaredCut)
  network <- abs(cor(data))^res$powerEstimate
  #save(network,file='result_wgcnaSoftThreshold.rda')
  
  network <- TOMsimilarity(network)
  save(network,file=paste0(outputpath,'result_wgcnaTOM.rda'))
}