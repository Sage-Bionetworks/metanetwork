#' Implements mrnet 
#' 
#' This function infers the interaction network of `data` using the MRNET algorithm.
#' 
#' @param data Required. A gene expression matrix with rows as sample IDs and
#' columns as Gene or feature IDs.
#' @param path Optional. String containing the path to the aracne compiled 
#' executable. (Default = NULL)
#' @param pval Optional. Cutoff p-value to determine a coexpressed edge. If one is
#' specified aracne will produce the file `aracneNetwork.csv` if less than 1 it 
#' will produce the file `aracneThresholdNetwork.csv`. (Default = 1)
#' @param outputpath Required. The output path to save the resulting coexpression
#' network
#'  
#' @return NULL. Saves a sparrow network object to paste0(`outputpath`,
#' `regressionFunction`,'mrnetNetwork.csv')
#'
#' @importFrom magrittr %>%
#' @export
mrnetWrapper = function(data,path=NULL,pval=1,outputpath){
  #library(parmigene)
  metanetwork::aracne(data,path,pval,outputpath)
  #library(data.table)
  #library(dplyr)
  if(pval==1){
    fileName <- paste0(outputpath,'aracneNetwork.csv')
  }else{
    fileName <- paste0(outputpath,'aracneThresholdNetwork.csv')
  }
  cat('fileName:',fileName,'\n')
  #load(fileName)
  #data.matrix(data.frame(data.table::fread('~/Desktop/sparrowZNetwork.csv',data.table=F),row.names=1))
  network <- data.table::fread(fileName,data.table=F) %>%
    data.frame(row.names=1) %>%
    data.matrix
  network <- network+t(network)
  gc()
  network <- parmigene::mrnet(data.matrix(network))
  #save(network,file=paste0(outputpath,'result_mrnet.rda'))
  network <- network*upper.tri(network)
  utils::write.csv(network,file=paste0(outputpath,'mrnetNetwork.csv'),quote=F)
}