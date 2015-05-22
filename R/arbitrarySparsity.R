arbitrarySparsity <- function(nedges,network){
  require(dplyr)
  #make the network into a matrix for sure (in case it is stored as a data frame or list)
  network <- network %>% as.matrix
  #set the diagonal to zero
  diag(network) <- 0
  #take the absolute value
  network <- network %>% abs
  #grab the upper diagonal
  networkVec <- network[network %>% upper.tri %>% which] %>% c
  #sort the upper diagonal, decreasing!
  sortedNetworkVec<-sort(networkVec,decreasing=T)
  #need to add a conditional if nedges exceeds the number of non-zero elements
  thres <- sortedNetworkVec[nedges]
  if(thres==0){
    return(NA)
  }else{
    return(network>thres)
  }
}