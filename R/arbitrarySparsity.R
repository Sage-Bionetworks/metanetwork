arbitrarySparsity <- function(network,nedges){
  diag(network) <- 0
  network <- abs(network)
  require(dplyr)
  networkVec <- c(network[network %>% upper.tri() %>% which])
  sortedNetworkVec<-sort(networkVec,decreasing=T)
  thres <- sortedNetworkVec[nedges]
  return(network>thres)
}