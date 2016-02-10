c3netWrapper  = function(data,path=NULL,pval=1,outputpath){
  library(c3net)
  network <- c3net(t(data))
  #save(network,file=paste0(outputpath,'result_c3net.rda'))
  write.csv(network,file=paste0(outputpath,'c3netNetwork.csv'),quote=F)
}