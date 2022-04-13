#' Ranks Consensus networks with ties.
#' 
#' This function ranks a list of consensus network objects and returns the best 
#' rank network. This version of rank consensus is compatible with ties.
#' 
#' @param networks Required. A list object containing an individual network as 
#' a list entry.
#'
#' @return The best rank consensus network.
#' @importFrom magrittr %>%
#' @export
rankConsensus2 <- function(networks){
  #library(bit64)
  #library(dplyr)
  
  aggregateRankFunction <- function(network, 
                                    upperTriIndices, 
                                    aggregateRank) {
    collapsedEdgeSet <- network[upperTriIndices]
    set.seed(123456)
    foo <- rank(-abs(collapsedEdgeSet), ties.method = "random") %>% 
      bit64::as.integer64()
    aggregateRank <- aggregateRank + foo
    return(aggregateRank)
  }
  
  upperTriIndices <- networks[[1]] %>% upper.tri %>% which
  cat("extracted upper triangular indices\n")
  
  lowerTriIndices <- networks[[1]] %>% lower.tri %>% which
  cat("extracted lower triangular indices\n")
  
  aggregateRank <- rep(0, length(upperTriIndices))
  
  for (i in 1:length(networks)) {
    cat("i:", i, "\n")
    aggregateRank <- aggregateRankFunction(networks[[i]], 
                                           upperTriIndices, 
                                           aggregateRank)
    gc()
  }
  
  cat("building final rank\n")
  print(aggregateRank[1:10])
  
  #library(bit64)
  print(aggregateRank[1:10])
  cat("make negative\n")
  aggregateRank <- -aggregateRank
  print(aggregateRank[1:10])
  cat("newway\n")
  finalRank <- bit64::rank.integer64(aggregateRank)
  cat("oldway\n")
  cat("renormalizing final rank\n")
  finalRank <- finalRank/max(finalRank)
  cat("turning into network\n")
  network <- matrix(0, nrow(networks[[1]]), ncol(networks[[1]]))
  colnames(network) <- colnames(networks[[1]])
  rownames(network) <- rownames(networks[[1]])
  network[upperTriIndices] <- finalRank
  network[lowerTriIndices] <- 0
  return(network)
}
