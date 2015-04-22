correlationBonferroni <- function(data,path=NULL){
  require(dplyr)
  thres <- 0.05/(ncol(data) %>% choose(2));
  pvals <- corPvalue(data)
  diag(pvals) <- 1
  #return(pvals<thres)
  network <- pvals<thres
  save(network,file='result_correlationBonferroni.rda')
}