# Function to find global modularity (Q)
compute.Modularity <- function(adj, mod, method = 'Newman1'){

  # Input
  #      adj = n x n adjacency matrix in the ltCMatrix format
  #      mod = n x 3 dimensional data frame with column names as Gene.ID, moduleNumber, and moduleLabel
  #      method = 'Newman1' or 'Newman2'
  
  # Output (list of following elements)
  #      Q = global modularity index
  
  # Error functions
  if(class(adj) != "ltCMatrix")
    stop('Adjacency matrix should be of class ltCMatrix')
  
  if(dim(adj)[1] != dim(adj)[2])
    stop('Adjacency matrix should be symmetric')
  
  if(dim(mod)[2] != 3)
    stop('Module label matrix should be a nx3 data frame')
  
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'upper', weighted = NULL, diag = F)
  adj = as.matrix(adj) + t(as.matrix(adj))
  
  # Get modules
  modules = mod$moduleNumber+1
  names(modules) = mod$Gene.ID
  
  if (method == 'Newman1'){
    Q = igraph::modularity(g, modules)
  } else {
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
    
    # Calculate global modularity
    Q = 0
    for (ci in names(comm)){
      Ein = edge.comm[ci,ci]
      Eout = sum(edge.comm[ci,], na.rm = T) - edge.comm[ci,ci]
      E = sum(edge.comm, na.rm = T)
      Q = Q + ((Ein/E)-((2*Ein + Eout)/(2*E))^2)
    }
    
  }
  return(Q)
}