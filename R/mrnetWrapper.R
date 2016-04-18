mrnetWrapper = function(data,path=NULL,pval=1,outputpath){
  library(parmigene)
  metanetwork::aracne(data,path,pval,outputpath)
  library(data.table)
  library(dplyr)
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
  write.csv(network,file=paste0(outputpath,'mrnetNetwork.csv'),quote=F)
}