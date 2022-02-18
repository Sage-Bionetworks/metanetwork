#' Find Global Modularity (Q)
#' 
#' This function finds Function to find global modularity (Q) of a upper triangular adjacency matrix.
#'
#' @param method Optional. If 'Newman1' is specified modularity is instituded in 
#' as-per the Clauset, Newman model in the function igraph::graph.modularity()
#' (Default = 'Newman1')
#'  
#' @inheritParams compute.LocalModularity
#'  
#' @return Q = modularity index.
#' 
#' @importFrom rlang .data
#' 
#' @export compute.Modularity
compute.Modularity <- function(adj, mod, method = 'Newman1'){
  # Error functions
  if(class(adj) != "matrix")
    stop('Adjacency matrix should be of class matrix')
  
  if(dim(adj)[1] != dim(adj)[2])
    stop('Adjacency matrix should be symmetric')
  
  if(!all(adj[lower.tri(adj)] == 0))
    stop('Adjacency matrix should be upper triangular')
  
  if(dim(mod)[2] != 3)
    stop('Module label matrix should be a nx3 data frame')
  
  # Get modules
  modules = mod$moduleNumber+1
  names(modules) = mod$Gene.ID
  
  if (method == 'Newman1'){
    # Convert to igraph graph object
    g = igraph::graph.adjacency(adj, mode = 'upper', weighted = TRUE, diag = F)
    Q = igraph::modularity(g, modules)
  } else {
    adj = adj + t(adj)
    
    # Get unique communities
    comm = plyr::dlply(mod, .variables = "moduleNumber", .fun = function(x){ unique(x$Gene.ID) }, .parallel = T)
    
    # Get number of edges between communities
    edge.comm = foreach::foreach(ci = names(comm), .packages = c('foreach', 'doParallel'), .combine = cbind) %dopar% {
      foreach::foreach(cj = names(comm), .combine = c) %dopar% {
        sum(adj[as.character(comm[[ci]]), as.character(comm[[cj]])], na.rm = T)
      }
    }
    rownames(edge.comm) = names(comm)
    colnames(edge.comm) = names(comm)
    
    # Calculate global modularity
    Q = foreach::foreach(ci=names(comm), .combine = c) %dopar% {
      Ein = edge.comm[ci,ci]
      Eout = sum(edge.comm[ci,], na.rm = T) - edge.comm[ci,ci]
      E = sum(edge.comm, na.rm = T)
      ((Ein/E)-((2*Ein + Eout)/(2*E))^2)
    }
    Q = sum(Q, na.rm = T)
  }
  return(Q)
}