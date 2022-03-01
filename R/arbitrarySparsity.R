#' This function Builds an arbitary sparse network
#' 
#' This introduces sparsity into a data frame or matrix with increasing frequency
#' as the `nedges` parameter is specified
#' 
#' @param nedges Required. Edges parameter to increase sparsity. 
#' @param network Required. A matrix or data frame of gene expression
#' 
#' @export 
#' @return a logical matrix of TRUE/FALSE where True indicates sparsity
#' 
arbitrarySparsity <- function(nedges,network){
  #require(dplyr)
  gc()
  cat('nedges: ',nedges,'\n')
  #make the network into a matrix for sure (in case it is stored as a data frame or list)
  #network <- network %>% as.matrix
  #set the diagonal to zero
  #diag(network) <- 0
  #take the absolute value
  network <- network %>% abs
  #grab the upper diagonal
  networkVec <- network[network %>% upper.tri %>% which] %>% c
  #sort the upper diagonal, decreasing!
  sortedNetworkVec<-sort(networkVec,decreasing=T)
  #need to add a conditional if nedges exceeds the number of non-zero elements
  thres <- sortedNetworkVec[nedges]
  gc()
  if(thres==0){
    return(NA)
  }else{
    return(network>=thres)
  }
}
