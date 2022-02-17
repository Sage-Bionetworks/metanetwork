#' Find Modules with Leading Eigen Edges
#' 
#' This function iteratively wraps a function that a tries to find densely 
#' connected subgraphs in a graph by calculating the leading non-negative 
#' eigenvector of the modularity matrix of the graph.
#' 
#' @inheritParams findModules.CFinder
#' 
#' @return  GeneModules = n x 3 dimensional data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
#' 
#' @importFrom magrittr %>%
#' @export
findModules.leading_eigen <- function(adj, nperm = 10, min.module.size = 30){
  # Error functions
  if(class(adj) != "matrix")
    stop('Adjacency matrix should be of class matrix')
  
  if(dim(adj)[1] != dim(adj)[2])
    stop('Adjacency matrix should be symmetric')
  
  if(!all(adj[lower.tri(adj)] == 0))
    stop('Adjacency matrix should be upper triangular')
  
  # Make adjacency matrix symmetric
  adj = adj + t(adj)
  adj[diag(adj)] = 0
  
  # Compute modules by permuting the labels nperm times
  all.modules = plyr::llply(1:nperm, .fun= function(i, adj, min.module.size){
    # Permute gene ordering
    ind = sample(1:dim(adj)[1], dim(adj)[1], replace = FALSE)
    adj1 = adj[ind,ind]
    
    # Find modules 
    mod = findModules.leading_eigen.once(adj1, min.module.size)
    
    # Compute local and global modularity
    adj1[lower.tri(adj1)] = 0
    Q = compute.Modularity(adj1, mod)
    Qds = compute.ModularityDensity(adj1, mod)
    
    return(list(mod = mod, Q = Q, Qds = Qds))
  }, adj, min.module.size)
  
  # Find the best module based on Q and Qds
  tmp = plyr::ldply(all.modules, function(x){
    data.frame(Q = x$Q, Qds = x$Qds)
  }) %>%
    dplyr::mutate(r = base::rank(Q)+base::rank(Qds))
  ind = which.max(tmp$r)
  
  mod = all.modules[[ind]]$mod
  
  return(mod)
}
