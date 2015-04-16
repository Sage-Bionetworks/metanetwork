correlationBonferroni <- function(data,path=NULL){
  thres <- 0.05/choose(ncol(data),2);
  pvals <- corPvalue(data)
  diag(pvals) <- 1
  return(pvals<thres)
}