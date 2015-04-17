arbitrarySparsity <- function(network,nedges){
  require(dplyr)
  diag(network) <- 0
  network <- network %>% abs
  require(dplyr)
  networkVec <- c(network[network %>% upper.tri() %>% which])
  sortedNetworkVec<-sort(networkVec,decreasing=T)
  thres <- sortedNetworkVec[nedges]
  return(network>thres)
}