deployNetworkSparsity <- function(network,sparsity){
  ###fxn to deploy a sparse network
  require(dplyr)
  return(sparsity %>% sapply(arbitrarySparsity,network))  
}