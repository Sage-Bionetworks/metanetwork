score.nodes <- function(g, G, h=3){
  
  # Get background genes from graph
  background.genes = igraph::V(g)$name
  
  # Calculate node degree for identifying global drivers
  node.degree = igraph::degree(g)
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = sd(node.degree, na.rm = T)
  
  # Find H-Layer neighborhood
  neighbor.nodes = lapply(1:h, function(hi,sg){
    igraph::ego(sg, order = hi, nodes = igraph::V(sg), mode = 'all', mindist = hi)
  }, g)
  
  # Find summed weight of current node based on h-layer neighbors
  node.scores = foreach::foreach (x = 1:length(neighbor.nodes)) %dopar% {
    score = sapply(neighbor.nodes[[x]], function(y, G){
      sum(G[names(y)] * (1/node.degree[names(y)]))
    }, G)
  } %>%
    do.call(cbind,.)
  node.scores = rowSums(node.scores * t(matrix(rep(1/c(1:h), length(G)), h, length(G))), na.rm = T)
  node.scores = node.scores + G
  names(node.scores) = background.genes
  
  return(node.scores)
}