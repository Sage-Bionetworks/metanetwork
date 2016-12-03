# Function to find modularity density (Qds)
compute.ModularityDensity <- function(adj, mod){
  
  # Input
  #      adj = n x n upper triangular adjacency in the matrix class format
  #      mod = n x 3 dimensional data frame with column names as Gene.ID, moduleNumber, and moduleLabel
  
  # Output
  #      Qds = module density
  
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
  adj = adj + t(adj)
  
  # Get unique communities
  comm = plyr::dlply(mod, .(moduleNumber), .fun = function(x){ unique(x$Gene.ID) }, .parallel = T)
  
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
  
  # Calculate local modularity
  Qds = foreach::foreach(ci=names(comm), .combine = c) %dopar% {
    if(edge.comm[ci,ci] != 0){
      Ein = edge.comm[ci,ci]
      Eout = sum(edge.comm[ci,], na.rm = T) - edge.comm[ci,ci]
      E = sum(edge.comm, na.rm = T)
      dc = 2*Ein/(mod.sz[ci] * (mod.sz[ci] - 1))
      
      Ecc.dcc = 0
      for (cj in names(comm)){
        if (ci != cj){
          Ecc.dcc = Ecc.dcc + edge.comm[ci, cj] ^2/(mod.sz[ci] * mod.sz[cj])
        }
      }
      Ecc.dcc = Ecc.dcc/ (2 * E)
      
      ((Ein/E) * dc) - ((2 * Ein + Eout) * dc/(2 * E))^2 - Ecc.dcc
    }
  }
  Qds = sum(Qds, na.rm = T)
  
  return(Qds)
}