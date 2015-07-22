#### Function to obtain basic network properties ####
networkProperties <- function(g,
                              metrics = c('density','diameter','avgNodeDegree','avgPathLength',
                                          'avgClusteringCoefficient','centralisation')){
  
  require(igraph)
  results <- list()
  
  # Density: Shows how sparse or dense a network is, for full network the value is 1 and for empty network the value is NaN(0/0)
  if (any(metrics %in% 'density')){  
    density <- graph.density(g)
    results <- c(results,list(density = density))
  }
  
  # Diameter: Maximum distance between all the vertex in network
  if (any(metrics %in% 'diameter')){  
    diameter <- diameter(g)
    results <- c(results,list(diameter = diameter))
  }
  
  # Average node degree: Average degree of all nodes in network
  if (any(metrics %in% 'avgNodeDegree')){  
    avgNodeDegree <- mean(degree(g))
    results <- c(results,list(avgNodeDegree = avgNodeDegree))
  }
  
  # Average path length: Average distance between any two vertex in network
  if (any(metrics %in% 'avgPathLength')){  
    avgPathLength <- average.path.length(g)
    results <- c(results,list(avgPathLength = avgPathLength))
  }
  
  # Average clustering coefficeint: Gives first hand information on network clusters
  if (any(metrics %in% 'avgClusteringCoefficient')){  
    V(g)$ClusteringCoefficeint <- transitivity(g,
                                               v=V(g),
                                               type = 'barrat')
    avgClusteringCoefficient <- mean(V(g)$ClusteringCoefficeint,na.rm=T)
    results <- c(results,list(avgClusteringCoefficient = avgClusteringCoefficient))
  }
  
  # Centralisation: Measures how well network has a star-like topology (closer to 1 more likely the topology is star-like)
  if (any(metrics %in% 'centralisation')){  
    centralisation <- vcount(g)/(vcount(g) - 2) * (max(degree(g))/(vcount(g) - 1) - graph.density(g))
    results <- c(results,list(centralisation = centralisation))
  }
  
  return(results)
}