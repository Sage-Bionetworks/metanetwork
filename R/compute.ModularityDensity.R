# Function to find modularity density (Qds)
compute.ModularityDensity <- function(adj, mod){
  
  # Input
  #      adj = n x n adjacency matrix in the ltCMatrix format
  #      mod = n x 3 dimensional data frame with column names as Gene.ID, moduleNumber, and moduleLabel
  
  # Output (list of following elements)
  #      Qds = modularity density
  
  # Error functions
  if(class(adj) != "ltCMatrix")
    stop('Adjacency matrix should be of class ltCMatrix')
  
  if(dim(adj)[1] != dim(adj)[2])
    stop('Adjacency matrix should be symmetric')
  
  if(dim(mod)[2] != 3)
    stop('Module label matrix should be a nx3 data frame')
  
  # Convert lsparseNetwork upper to symmetric
  adj = as.matrix(adj) + t(as.matrix(adj))
  
  # Get unique communities
  comm = plyr::dlply(mod, .(moduleNumber), .fun = function(x){ unique(x$Gene.ID) })
  
  # Get number of edges between communities
  edge.comm = matrix(0, length(comm), length(comm))
  rownames(edge.comm) = names(comm)
  colnames(edge.comm) = names(comm)
  
  for(ci in names(comm)){
    for(cj in names(comm)){
      edge.comm[ci, cj] = sum(adj[as.character(comm[[ci]]), as.character(comm[[cj]])], na.rm = T)
    }
  }
  
  # Get size of each modules
  mod.sz = sapply(comm, length)
  
  # Calculate local modularity
  Qds = 0
  for (ci in names(comm)){
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
      
      Qds = Qds + ((Ein/E) * dc) - ((2 * Ein + Eout) * dc/(2 * E))^2 - Ecc.dcc
    }
  }
  
  return(Qds)
}