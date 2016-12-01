# Function to find local modularity (NQ)
compute.LocalModularity <- function(adj, mod){
  
  # Input
  #      adj = n x n adjacency matrix in the ltCMatrix format
  #      mod = n x 3 dimensional data frame with column names as Gene.ID, moduleNumber, and moduleLabel
  
  # Output (list of following elements)
  #      NQ = local modularity index
  
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
    
  # Calculate local modularity
  NQ = 0
  for (ci in names(comm)){
    if(edge.comm[ci,ci] != 0){
      Ein = edge.comm[ci,ci]
      Eout = sum(edge.comm[ci,], na.rm = T) - edge.comm[ci,ci]
      ind = which(edge.comm[ci,] != 0)
      Eneighbor = sum(edge.comm[ind,ind], na.rm = T) 
      NQ = NQ + ((Ein/Eneighbor)-((2*Ein + Eout)/(2*Eneighbor))^2)
    }
  }
  
  return(NQ)
}