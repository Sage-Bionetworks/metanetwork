wgcnaSoftThreshold <- function(data,path=NULL,outputpath){
  library(WGCNA)
  res <- WGCNA::pickSoftThreshold(data,RsquaredCut=0.80)
  network <- abs(cor(data))^res$powerEstimate
  #save(network,file='result_wgcnaSoftThreshold.rda')
  save(network,file=paste0(outputpath,'result_wgcnaSoftThreshold.rda'))
}