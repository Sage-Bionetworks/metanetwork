rankConsensus <- function(networks){
  library(bit64)
  library(dplyr)
  aggregateRankFunction <- function(network,upperTriIndices,aggregateRank){
    collapsedEdgeSet <- network[upperTriIndices]
    foo <- rank(-abs(collapsedEdgeSet),ties.method='min') %>% as.integer64
    #foo <- rank(-abs(collapsedEdgeSet),ties.method='min') %>% as.integer
    aggregateRank <- aggregateRank + foo
    return(aggregateRank)
  }
  upperTriIndices <- networks[[1]] %>%
                     upper.tri %>%
                     which
  cat('extracted upper triangular indices\n')
  lowerTriIndices <- networks[[1]] %>%
                     lower.tri %>%
                     which
  cat('extracted lower triangular indices\n')
  aggregateRank <- rep(0,length(upperTriIndices))
  for (i in 1:length(networks)){
    cat('i:',i,'\n')
    aggregateRank <- aggregateRankFunction(networks[[i]],upperTriIndices,aggregateRank)
    gc()
  }
  save(aggregateRank,file='/shared/CRANIO/aggregateRank.rda')
  cat('building final rank\n')
  finalRank <- rank(-aggregateRank,ties.method = 'min')
  cat('renormalizing final rank\n')
  finalRank <- finalRank/max(finalRank)
  cat('turning into network\n')
  network <- matrix(0,nrow(networks[[1]]),ncol(networks[[1]]))
  colnames(network) <- colnames(networks[[1]])
  rownames(network) <- rownames(networks[[1]])
  network[upperTriIndices] <- finalRank
  network[lowerTriIndices] <- 0
  return(network)
}
