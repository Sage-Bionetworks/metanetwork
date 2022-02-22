#' This function deploys arbitary sparsity on a list of networks
#' 
#' This introduces sparsity into a list of data frames or matrices with the 
#' `arbitrarySparsity()` function with increasing frequency
#' as the `sparsity` parameter is specified
#' 
#' @param sparsity Required. Edges parameter to increase sparsity. 
#' @param network Required. A matrix or data frame of gene expression
#' 
#' @export 
#' @return a list of logical matrices of TRUE/FALSE where True indicates sparsity
#' 
deployNetworkSparsity <- function(network,sparsity){
  ###fxn to deploy a sparse network
  #require(dplyr)
  return(sparsity %>% lapply(arbitrarySparsity,network))
}