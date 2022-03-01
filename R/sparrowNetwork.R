#' Runs wraps vbsr across a matrix 
#' 
#' This function wraps variable bays spike regression calling all vsbr on all genes 
#' pairwise across an expression matrix
#' 
#' @param data Required. A .
#' 
#' @return A coexpression matrix
#' @export
sparrowNetwork <- function(data){
  #library(vbsr)
  thr <- stats::qchisq(0.05/choose(ncol(data),2),1,lower.tail=F)
  network <- matrix(0,ncol(data),ncol(data))
  rownames(network) <- colnames(data)
  colnames(network) <- colnames(data)
  for (i in 1:ncol(data)){
    network[i,-i] <- vbsr::vbsr(data[,i],data[,-i])$z
  }
  network <- network/2+t(network/2)
  #network <- network^2>thr
  return(network)
}