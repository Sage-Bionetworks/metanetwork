#' Runs WGCNA TOMsimilarity
#' 
#' Runs WGCNA::TOMsimilarity over a gene expression matrix.
#' 
#' @inheritParams wgcnaSoftThreshold
#' 
#' @return NULL. Writes coexpression network named wgcnaTopologicalOverlapMatrixNetwork.csv
#' to `outputpath`
#' 
#' @export
wgcnaTOM <- function(data,path=NULL,pval=1,outputpath,RsquaredCut=.80,defaultNaPower=6){
  #library(WGCNA)
  res <- WGCNA::pickSoftThreshold(data,RsquaredCut=RsquaredCut)
  if(is.na(res$powerEstimate)){
    res$powerEstimate<-defaultNaPower
  }
  network <- abs(stats::cor(data))^res$powerEstimate
  utils::write.csv(network*upper.tri(network),file=paste0(outputpath,'wgcnaSoftThresholdNetwork.csv')) 
  gc()
  cat(res$powerEstimate,'\n',sep='',file=paste0(outputpath,'wgcnaPowerEstimate.txt'))
  #save(network,file='result_wgcnaSoftThreshold.rda')
  print(res$powerEstimate)
  cn <- colnames(network)
  network <- WGCNA::TOMsimilarity(network)
  colnames(network) <- cn
  rownames(network) <- cn
  #save(network,file=paste0(outputpath,'result_wgcnaTOM.rda'))
  network <- network*upper.tri(network)
  utils::write.csv(network,file=paste0(outputpath,'wgcnaTopologicalOverlapMatrixNetwork.csv'),quote=F)
}
