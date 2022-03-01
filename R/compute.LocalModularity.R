#' Find Local Modularity (NQ)
#' 
#' This function finds local modularity (NQ) in an upper triangular adjacency matrix.
#'
#' @param adj Required. An n x n upper triangular adjacency in the matrix class 
#' format.
#' @param mod Required. An n x 3 dimensional data frame with column names as 
#' Gene.ID, moduleNumber, and moduleLabel.
#'  
#' @return NQ = local modularity index.
#' 
#' @export compute.LocalModularity
#' 
compute.LocalModularity <- function(adj, mod){
  # Error functions
  if(class(adj) != "matrix")
    stop('Adjacency matrix should be of class matrix')
  
  if(dim(adj)[1] != dim(adj)[2])
    stop('Adjacency matrix should be symmetric')
  
  if(!all(adj[lower.tri(adj)] == 0))
    stop('Adjacency matrix should be upper triangular')
  
  if(dim(mod)[2] != 3)
    stop('Module label matrix should be a nx3 data frame')
  
  # Convert lsparseNetwork upper adj matrix to graph
  g = igraph::graph.adjacency(adj, mode = 'upper', weighted = TRUE, diag = F)
  rownames(mod) = mod$Gene.ID
  igraph::V(g)$moduleNumber = paste0('mod.',mod[igraph::V(g)$name, 'moduleNumber'])
  
  rm(list = c('adj', 'mod'))
  gc()
  
  # Get number of edges between communities
  edge.comm = foreach::foreach(ci=unique(igraph::V(g)$moduleNumber), 
                               .packages = c('foreach', 'doParallel'), 
                               .combine = cbind, 
                               .export = c('g')) %dopar% {
                                 foreach::foreach(cj=unique(igraph::V(g)$moduleNumber), 
                                                  .combine = c,
                                                  .packages = c('foreach', 'doParallel', 'dplyr'),
                                                  .export = c('g', 'ci')) %dopar% {
                                                    gi = igraph::induced_subgraph(g, vids = which(igraph::V(g)$moduleNumber == ci))
                                                    gj = igraph::induced_subgraph(g, vids = which(igraph::V(g)$moduleNumber == cj))
                                                    igraph::ecount(igraph::intersection(gi,gj))
      }
  }
  edge.comm = data.frame(edge.comm)
  rownames(edge.comm) = unique(igraph::V(g)$moduleNumber)
  colnames(edge.comm) = unique(igraph::V(g)$moduleNumber)
  
  rm(list = c('g'))
  gc()
  
  # Calculate local modularity
  NQ = foreach::foreach(ci=rownames(edge.comm), 
                        .combine = c,
                        .export = c('edge.comm')) %dopar% {
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