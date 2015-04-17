arbitrarySparsity <- function(network,nedges){
  require(dplyr)
  network <- network %>% as.matrix
  diag(network) <- 0
  network <- network %>% abs
  require(dplyr)
  networkVec <- network[network %>% upper.tri %>% which] %>% c
  sortedNetworkVec<-sort(networkVec,decreasing=T)
  thres <- sortedNetworkVec[nedges]
  return(network>thres)
}