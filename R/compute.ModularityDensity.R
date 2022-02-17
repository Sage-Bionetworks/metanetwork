#' Find Global Modularity (Qds)
#' 
#' This function finds calculates modularity density.
#'
#' @inheritParams compute.LocalModularity
#'  
#' @return Qds = module density.
#' 
#' @export
compute.ModularityDensity <- function(adj, mod){
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
  edge.comm = foreach::foreach(ci = unique(igraph::V(g)$moduleNumber), 
                               .packages = c('foreach', 'doParallel'), 
                               .combine = cbind, 
                               .export = c('g')) %dopar% {
                                 foreach::foreach(cj = unique(igraph::V(g)$moduleNumber), 
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
  
  # Get size of each modules
  mod.sz = summary(factor(igraph::V(g)$moduleNumber))
  
  rm(list = c('g'))
  gc()
  
  # Calculate local modularity
  Qds = foreach::foreach(ci = rownames(edge.comm),
                         .combine = c,
                         .export = c('edge.comm', 'mod.sz')) %dopar% {
                           if(edge.comm[ci,ci] != 0){
                             Ein = edge.comm[ci,ci]
                             Eout = sum(edge.comm[ci,], na.rm = T) - edge.comm[ci,ci]
                             E = sum(edge.comm, na.rm = T)
                             dc = 2*Ein/(mod.sz[ci] * (mod.sz[ci] - 1))
                             
                             Ecc.dcc = 0
                             
                             for (cj in rownames(edge.comm)){
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