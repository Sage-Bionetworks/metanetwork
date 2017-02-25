# Function to find qulaity metrics of all modules
compute.ModuleQualityMetric <- function(adj, mod){
  
  # Input
  #      adj = n x n upper triangular adjacency in the matrix class format
  #      mod = n x 3 dimensional data frame with column names as Gene.ID, moduleNumber, and moduleLabel
  
  # Output (list of following elements)
  #      metrics = data frame of module quality metrics
  
  # Error functions
  if(class(adj) != "matrix")
    stop('Adjacency matrix should be of class matrix')
  
  if(dim(adj)[1] != dim(adj)[2])
    stop('Adjacency matrix should be symmetric')
  
  if(!all(adj[lower.tri(adj)] == 0))
    stop('Adjacency matrix should be upper triangular')
  
  if(dim(mod)[2] != 3)
    stop('Module label matrix should be a nx3 data frame')
  
  # Convert upper to symmetric matrix
  adj = adj + t(adj)
  
  # Get unique communities
  comm = plyr::dlply(mod, .variables = 'moduleNumber', .fun = function(x){ unique(x$Gene.ID) }, .parallel = T)
  
  
  # Get number of edges between communities
  edge.comm = foreach::foreach(ci = names(comm), .packages = c('foreach', 'doParallel'), .combine = cbind) %dopar% {
    foreach::foreach(cj = names(comm), .combine = c) %dopar% {
      sum(adj[as.character(comm[[ci]]), as.character(comm[[cj]])], na.rm = T)
    }
  }
  rownames(edge.comm) = names(comm)
  colnames(edge.comm) = names(comm)
  
  # Get size of each modules
  mod.sz = sapply(comm, length)
  
  # Intra edges
  IE = data.frame(moduleNumber = names(comm),
                  IE = diag(edge.comm))
  
  # Intra density
  ID = data.frame(moduleNumber = names(comm),
                  ID = 2*diag(edge.comm)/(mod.sz * (mod.sz - 1)))
  
  # Contraction
  CNT = data.frame(moduleNumber = names(comm),
                   CNT = 2*diag(edge.comm)/mod.sz)
  
  # Boundary Edges
  BE = data.frame(moduleNumber = names(comm),
                  BE = rowSums(edge.comm, na.rm = T) - diag(edge.comm))
  
  # Expansion
  EXP = data.frame(moduleNumber = names(comm),
                   EXP = BE$BE/mod.sz)
  
  # Conductance
  CND = data.frame(moduleNumber = names(comm),
                   CND = BE$BE/(2*IE$IE + BE$BE))
  
  # Fitness function
  FIT = data.frame(moduleNumber = names(comm),
                   FIT = IE$IE/(IE$IE + BE$BE))
  
  # Average modularity degree
  AMD = data.frame(moduleNumber = names(comm),
                   AMD = (2*IE$IE - BE$BE)/mod.sz)
  
  metrics = plyr::join_all(list(IE, ID, CNT, BE, EXP, CND, FIT, AMD)) %>%
    dplyr::inner_join(mod %>% 
                        dplyr::select(moduleNumber, moduleLabel) %>%
                        dplyr::mutate(moduleNumber = factor(moduleNumber)) %>%
                        unique, .)
  return(metrics)
}