#### Function to measure different node ranking metrics ####
nodeRankings <- function(g, # graph for which network properties are to be calculated
                         metrics = 'totalDegree' # character vector of network properties, 
# Other options are:
    # inDegree
    # outDegree
    # totalDegree
    # betwenness
    # clustCoefficient
    # closeness
    # eigen
    # pageRank
    # nearNeighbor
                         ){
  
  require(igraph)
  results <- data.frame(row.names=V(g)$name)
  
  # In degree
  if (any(metrics %in% 'inDegree')){
    print('Calculating in degree...')
    results$inDegree <- degree(g,
                               v=V(g), 
                               mode = 'in', 
                               loop = T,
                               normalized = T)
  }
  
  # Out degree
  if (any(metrics %in% 'outDegree')){
    print('Calculating out degree...')
    results$outDegree <- degree(g,
                                v=V(g), 
                                mode = 'out', 
                                loop = T,
                                normalized = T)
  }
  
  # Total degree: An important node is involved in a large number of interactions
  if (any(metrics %in% 'totalDegree')){
    print('Calculating total degree...')
    results$totalDegree <- degree(g,
                                  v=V(g), 
                                  mode = 'all', 
                                  loop = T,
                                  normalized = T)
  }
  
  # Betwenness centrality: An important node will lie on a high proportion of paths between other nodes in the network.
  if (any(metrics %in% 'betwenness')){
    print('Calculating betweenness...')
    results$betweenness <- betweenness(g,
                                       v= V(g),
                                       normalized = T)
  }
  
  # Clustering coefficeint
  if (any(metrics %in% 'clustCoefficient')){
    print('Calculating clustering coefficient...')
    results$clustCoefficient <- transitivity(g,
                                             v=V(g),
                                             type = 'barrat')
  }
  
  # Nearest Neighbor degree
  if (any(metrics %in% 'nearNeighbor')){
    print('Calculating k nearest neighborhood degree...')
    results$nearNeighbor <- graph.knn(g,
                                      vids=V(g))$knn
  }
  
  # Page Rank:  An important node is likely to receive more connections from other nodes.
  if (any(metrics %in% 'pageRank')){
    print('Calculating in page rank...')
    results$pageRank <- page.rank(g,
                                  vids = V(g))$vector
  }
  
  # Eigen centrality: An important node is connected to important neighbours.
  if (any(metrics %in% 'eigen')){
    print('Calculating eigen centrality...')
    results$eigen <- alpha.centrality(g,
                                      nodes=V(g))
  }
  
  # Closeness centrality: An important node is typically “close” to, and can communicate quickly with, the other nodes in the network
  if (any(metrics %in% 'closeness')){
    print('Calculating closeness centrality ...')
    results$closeness <- closeness(g,
                                   v= V(g),
                                   mode = 'all',
                                   normalized = T)
  }
  
  return(results)
}