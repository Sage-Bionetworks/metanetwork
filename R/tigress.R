#' Runs tigress on an expression matrix
#' 
#' This function implements Trustful Inference of Gene REgulation with Stability 
#' Selection (TIGRESS) algoritm..
#' 
#' @param y Required. A vector of gene expression values.
#' @param x Required. A gene expression matrix.
#'  
#' @return Vector of coexpression values of gene Y to columns of X
#' @export
tigress<- function(y,x){
  alpha = 0.2;
  L = 5;
  R = 100;
  n <- length(y);
  #require(lars);
  #indexMat <- matrix(0,R,L)
  indexMat <- matrix(0,L,ncol(x))
  for(i in 1:floor(R/2)){
    indexVec <- sample(1:n,n);
    xr1 <- t( t( x[ indexVec[1:floor(n/2)],]) * stats::runif(ncol(x),alpha,1) );
    xr2 <- t(t(x[indexVec[(floor(n/2)+1):n],])*stats::runif(ncol(x),alpha,1));
    result1 <- lars::lars(x=xr1,y=y[indexVec[1:floor(n/2)]],type='lar',max.steps=L,use.Gram=FALSE)
    w1<-rev(order(result1$entry,decreasing=T)[1:L])
    indexMat[,w1][lower.tri(indexMat[,w1],diag=T)] <- indexMat[,w1][lower.tri(indexMat[,w1],diag=T)] + 1;
    result2 <- lars::lars(x=xr2,y=y[indexVec[(floor(n/2)+1):n]],type='lar',max.steps=L,use.Gram=FALSE)
    w2<-rev(order(result2$entry,decreasing=T)[1:L])
    indexMat[,w2][lower.tri(indexMat[,w2],diag=T)] <- indexMat[,w2][lower.tri(indexMat[,w2],diag=T)] + 1;
  }
  return(colMeans(indexMat/R))
}

