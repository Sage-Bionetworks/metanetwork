correlationBonferroni <- function(data,path=NULL,outputpath){
  require(dplyr)
  thres <- 0.05/(ncol(data) %>% choose(2));
  pvals <- corPvalue(data)
  diag(pvals) <- 1
  #return(pvals<thres)
  network <- pvals<thres
  cat(paste('correlationBonferroni',sum(network)/2,sep=','),'\n',file=paste0(outputpath,'sparsity.csv'),sep='',append=TRUE)
  save(network,file=paste0(outputpath,'result_correlationBonferroni.rda'))
}