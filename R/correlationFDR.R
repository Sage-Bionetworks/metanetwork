correlationFDR <- function(data,path=NULL,fdr=0.05,outputpath){
  #thres <- 0.05/choose(ncol(data),2);
  require(dplyr)
  pvals <- data %>% corPvalue
  diag(pvals) <- 1
  thres <- pvals[pvals %>% upper.tri %>% which] %>% fdrThres(fdr=fdr)
  #return(pvals<thres)
  network <- pvals < thres
  cat(paste('correlationFDR',sum(network)/2,sep=','),'\n',file=paste0(outputpath,'sparsity.csv'),sep='',append=TRUE)
  save(network,file=paste0(outputpath,'result_correlationFDR.rda'))
}