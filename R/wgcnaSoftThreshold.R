#' Runs WGCNA
#' 
#' Runs WGCNA::pickSoftThreshold over a gene expression matrix.
#' 
#' @param data Required. Gene expression data in a matrix or data frame. Rows 
#' correspond to samples and columns to genes.
#' @param outputpath Required. The file path to write the resulting coexpression 
#' network. 
#' @param path Optional. Not utilized in this function. (Default = NULL)
#' @param pval Optional. Not utilized in this function. (Default = 1)
#' @param RsquaredCut Optional. Desired minimum scale free topology fitting index 
#' R^2. (Default = 0.80)
#' @param defaultNaPower Optional. The power to rais the abs(cor(data)) matrix to
#' if thee power estimate from WGCNA::pickSoftThreshold() is too low (Default = 6).
#'  
#' @return NULL. Writes coexpression network named wgcnaSoftThresholdNetwork.csv
#' to `outputpath`
#' 
#' @export
wgcnaSoftThreshold <- function(data,path=NULL,pval=1,outputpath,RsquaredCut=.80,defaultNaPower=6){
  #library(WGCNA)
  res <- WGCNA::pickSoftThreshold(data,RsquaredCut=RsquaredCut)
  if(is.na(res$powerEstimate)){
    res$powerEstimate<-defaultNaPower
  }
  
  network <- abs(cor(data))^res$powerEstimate

  #save(network,file=paste0(outputpath,'result_wgcnaST.rda'))
  network <- network*upper.tri(network)
  write.csv(network,file=paste0(outputpath,'wgcnaSoftThresholdNetwork.csv'),quote=F)
}
