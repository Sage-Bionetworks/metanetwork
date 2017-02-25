# Function to find local modularity (NQ)
compute.LocalModularity <- function(adj, mod){
  
  # Input
  #      adj = n x n upper triangular adjacency in the matrix class format
  #      mod = n x 3 dimensional data frame with column names as Gene.ID, moduleNumber, and moduleLabel
  
  # Output (list of following elements)
  #      NQ = local modularity index
  
  # Error functions
  if(class(adj) != "matrix")
    stop('Adjacency matrix should be of class matrix')
  
  if(dim(adj)[1] != dim(adj)[2])
    stop('Adjacency matrix should be symmetric')
  
  if(!all(adj[lower.tri(adj)] == 0))
    stop('Adjacency matrix should be upper triangular')
  
  if(dim(mod)[2] != 3)
    stop('Module label matrix should be a nx3 data frame')
  
  # Convert lsparseNetwork upper to symmetric
  adj = as.matrix(adj) + t(as.matrix(adj))
  
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
  
  # Calculate local modularity
  NQ = foreach::foreach(ci = names(comm), .combine = c) %dopar% {
    if(edge.comm[ci,ci] != 0){
      Ein = edge.comm[ci,ci]
      Eout = sum(edge.comm[ci,], na.rm = T) - edge.comm[ci,ci]
      ind = which(edge.comm[ci,] != 0)
      Eneighbor = sum(edge.comm[ind,ind], na.rm = T) 
      ((Ein/Eneighbor)-((2*Ein + Eout)/(2*Eneighbor))^2)
    }
  }
  NQ = sum(NQ)
  
  return(NQ)
}