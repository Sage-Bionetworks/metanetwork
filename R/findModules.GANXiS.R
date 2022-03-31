#' Find Modules with GANXiS
#' 
#' This function to get modules from network adjacency matrix using GANXiS 
#' community detection algorithm v3.0.2.
#' 
#' @inheritParams findModules.CFinder
#' 
#' @return  GeneModules = n x 3 dimensional data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
#' 
#' @importFrom magrittr %>%
#' @export
findModules.GANXiS <- function(adj, path, nperm = 10, min.module.size = 30){
  
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
  all.modules = plyr::llply(1:nperm, .fun= function(i, adj, path, min.module.size){
    # Permute gene ordering
    ind = sample(1:dim(adj)[1], dim(adj)[1], replace = FALSE)
    adj1 = adj[ind,ind]
    
    mod = NA
    Q = NA
    Qds = NA
    tryCatch({
      # Find modules 
      mod = findModules.GANXiS.once(adj1, path, min.module.size)
    
      # Compute local and global modularity
      adj1[lower.tri(adj1)] = 0
      Q = compute.Modularity(adj1, mod)
      Qds = compute.ModularityDensity(adj1, mod)
    }, error = function(e){
      mod = NA; Q = NA; Qds = NA;
    })
    
    return(list(mod = mod, Q = Q, Qds = Qds))
  }, adj, path, min.module.size)
  
  # Find the best module based on Q and Qds
  tmp = plyr::ldply(all.modules, function(x){
    data.frame(Q = x$Q, Qds = x$Qds)
  }) %>% stats::na.omit %>% dplyr::mutate(r = base::rank(Q, na.last = FALSE)+base::rank(Qds, na.last = FALSE))
  ind = which.max(tmp$r)
  
  mod = all.modules[[ind]]$mod
  
  return(mod)
}
