rankConsensus <- function(networks){
  #library(bit64)
  library(dplyr)
  aggregateRankFunction <- function(network,upperTriIndices,aggregateRank){
    collapsedEdgeSet <- network[upperTriIndices]
    #foo <- rank(-abs(collapsedEdgeSet),ties.method='min') %>% as.integer64
    foo <- rank(-abs(collapsedEdgeSet),ties.method='min') %>% as.integer
    aggregateRank <- aggregateRank + foo
    return(aggregateRank)
  }
  upperTriIndices <- networks[[1]] %>%
                     upper.tri %>%
                     which

  lowerTriIndices <- networks[[1]] %>%
                     lower.tri %>%
                     which

  aggregateRank <- rep(0,length(upperTriIndices))
  for (i in 1:length(networks)){
    aggregateRank <- aggregateRankFunction(networks[[i]],upperTriIndices,aggregateRank)
    gc()
  }

  finalRank <- rank(-aggregateRank,ties.method = 'min')
  finalRank <- finalRank/max(finalRank)
  network <- matrix(0,nrow(networks[[1]]),ncol(networks[[1]]))
  colnames(network) <- colnames(networks[[1]])
  rownames(network) <- rownames(networks[[1]])
  network[upperTriIndices] <- finalRank
  network[lowerTriIndices] <- 0
  return(network)
}
