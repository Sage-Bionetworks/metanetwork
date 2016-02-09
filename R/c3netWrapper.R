c3netWrapper  = function(data,path=NULL,pval=1,outputpath){
  library(c3net)
  network <- c3net(t(data))
  save(network,file='result_c3net.rda')
}