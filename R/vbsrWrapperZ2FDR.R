#vbsrWrapper
vbsrWrapperZ2FDR <- function(y,x,fdr=0.05,...){
  require(vbsr)
  result <- vbsr(y=y,X=x,...)$z
  pval <- pchisq(result^2,1,lower.tail=F)  
  if(is.null(FDR)){
    thres <- 0.05/ncol(x);
  }else{
    thres <- fdrThres(pval,fdr = fdr)
  }
  if(sum(pval<thres)>0){
    newz <- fastlm(y,x[,pval<thres])
    result[pval<thres] <- newz
  }
  return(result)
}
