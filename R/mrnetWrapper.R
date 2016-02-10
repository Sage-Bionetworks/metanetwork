mrnetWrapper = function(data,path=NULL,pval=1,outputpath){
  library(parmigene)
  metanetwork::aracne(data,path,pval,outputpath)
  library(data.table)
  if(pval==1){
    fileName <- paste0(outputpath,'aracneNetwork.csv')
  }else{
    fileName <- paste0(outputpath,'aracneThresholdNetwork.csv')
  }
  cat('fileName:',fileName,'\n')
  #load(fileName)
  data.table::fread(fileName,data.table=F)
  network <- parmigene::mrnet(network)
  #save(network,file=paste0(outputpath,'result_mrnet.rda'))
  write.csv(network,file=paste0(outputpath,'mrnetNetwork.csv'),quote=F)
}