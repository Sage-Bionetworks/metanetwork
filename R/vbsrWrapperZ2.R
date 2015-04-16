#vbsrWrapper
vbsrWrapperZ2 <- function(y,x,...){
  require(vbsr)
  result <- vbsr(y=y,X=x,...)$z
  thres <- 0.05/ncol(x);
  pval <- pchisq(result^2,1,lower.tail=F)
  if(sum(pval<thres)>0){
    newz <- fastlm(y,x[,pval<thres])
    result[pval<thres] <- newz
  }
  return(result)
}
