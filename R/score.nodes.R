# Function to score network nodes based on its neighborhood scores
score.nodes <- function(g, G, h=3, mode = 'all'){
  
  # Input
  #      g = an igraph object with n vertices
  #      G = a named vector of node scores
  #      h = neighborhood search distance (h nodes away from current node)
  #      mode = one of c("all", "out", "in", "total"). Character string, “out” for out-degree, “in” for in-degree or “all” for the sum of the two. For undirected graphs this argument is ignored.
  
  # Output
  #      node.scores = n x 1 dimensional vector of node scores based on its neighborhood
  
  # Error checking
  if(!igraph::is.igraph(g))
    stop('Input has to be an igraph object')
  
  if(!is.numeric(G))
    stop('G has to be a named numeric vector')
  
  if(length(intersect(igraph::V(g)$name,names(G))) == 0)
    stop('Names of G has to match adjacency names')
  
  # Get background genes from graph
  background.genes = igraph::V(g)$name
  
  # Assign zero scores to nodes without any value
  G = G[intersect(background.genes, names(G))]
  G.not = rep(0, length(setdiff(background.genes, names(G))))
  names(G.not) = setdiff(background.genes, names(G))
  G = c(G, G.not)
  
  # Calculate node degree for identifying global drivers
  node.degree = igraph::strength(g, mode = mode, loops = FALSE)
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = sd(node.degree, na.rm = T)
  
  # Find H-Layer neighborhood
  neighbor.nodes = lapply(1:h, function(hi,sg){
    igraph::ego(sg, order = hi, nodes = igraph::V(sg), mode = mode, mindist = hi)
  }, g)
  
  # Find summed weight of current node based on h-layer neighbors
  node.scores = foreach::foreach (x = 1:length(neighbor.nodes), .combine = cbind) %dopar% {
    score = sapply(neighbor.nodes[[x]], function(y, G){
      sum(G[names(y)] * (1/node.degree[names(y)]))
    }, G)
  } 
  node.scores = rowSums(node.scores * t(matrix(rep(1/c(1:h), length(G)), h, length(G))), na.rm = T)
  node.scores = node.scores + G
  names(node.scores) = background.genes
  
  return(node.scores)
}