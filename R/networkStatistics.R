
networkStatistics <- function(network){
  #function to extract network statistics from data
  #degree
  #stress
  #betweeness centrality
  #topological coefficient
  #shared neighbors
  #clustering coefficient
  #closeness centrality
  netStat <- list()
  netStat$hubs <- rowSums(network)
  return(netStat)  
}
