#' Runs vbsr with a 2Z cutoff for a gene across a matrix 
#' 
#' This function wraps variable bays spike regression of a genes expression across
#' a matrix of genes expressed in the same samples with a 2Z cutoff value.
#' 
#' @param y Required. response variable. Normally distributed errors for 
#' family="normal". For family="binomial" should be coded as a vector of 0's and 1's.
#' @param x Required. Design matrix, an n x m matrix, with rows as observations.
#' @param fdr Optional. FDR threshold cut off for edge determination. NULL results
#' in a cutoff of 0.05 (Default = NULL)
#' 
#' @return A coexpression value
#' @export
sparrow2Z <- function(y,x,fdr=NULL,...){
  #require(vbsr)
  result <- vbsr::vbsr(y=y,X=x,...)$z
  pval <- stats::pchisq(result^2,1,lower.tail=F)  
  if(is.null(fdr)){
    thres <- 0.05/ncol(x);
  }else{
    thres <- utilityFunctions::fdrThres(pval,fdr = fdr)
  }
  if(sum(pval<thres)>0){
    newz <- utilityFunctions::fastlm(y,x[,pval<thres])
    result[pval<thres] <- newz
  }
  return(result)
}
