mrnetWrapper = function(data,path=NULL,pval=1,outputpath){
  library(parmigene)
  metanetwork::aracne(data,path,pval,outputpath)
  if(pval==1){
    fileName <- paste0(outputpath,'result_aracneFull.rda')
  }else{
    fileName <- paste0(outputpath,'result_aracne.rda')
  }
  load(fileName)
  network <- parmigene::mrnet(network)
  save(network,file='result_mrnet.rda')
}