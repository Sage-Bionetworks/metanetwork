#' Find  Modularity Quality
#' 
#' This function finds qulaity metrics of all modules
#'
#' @inheritParams compute.LocalModularity
#'  
#' @return metrics = data frame of module quality metrics.
#' 
#' @export
compute.ModuleQualityMetric <- function(adj, mod){
  
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
  edge.comm = data.matrix(edge.comm)
  
  # Get size of each modules
  mod.sz = summary(factor(igraph::V(g)$moduleNumber))
  
  rm(list = c('g'))
  gc()
  
  # Intra edges
  IE = data.frame(moduleNumber = rownames(edge.comm), IE = diag(edge.comm))
  
  # Intra density
  ID = data.frame(moduleNumber = rownames(edge.comm),
                  ID = 2*diag(edge.comm)/(mod.sz * (mod.sz - 1)))
  
  # Contraction
  CNT = data.frame(moduleNumber = rownames(edge.comm),
                   CNT = 2*diag(edge.comm)/mod.sz)
  
  # Boundary Edges
  BE = data.frame(moduleNumber = rownames(edge.comm),
                  BE = rowSums(edge.comm, na.rm = T) - diag(edge.comm))
  
  # Expansion
  EXP = data.frame(moduleNumber = rownames(edge.comm),
                   EXP = BE$BE/mod.sz)
  
  # Conductance
  CND = data.frame(moduleNumber = rownames(edge.comm),
                   CND = BE$BE/(2*IE$IE + BE$BE))
  
  # Fitness function
  FIT = data.frame(moduleNumber = rownames(edge.comm),
                   FIT = IE$IE/(IE$IE + BE$BE))
  
  # Average modularity degree
  AMD = data.frame(moduleNumber = rownames(edge.comm),
                   AMD = (2*IE$IE - BE$BE)/mod.sz)
  
  metrics = plyr::join_all(list(IE, ID, CNT, BE, EXP, CND, FIT, AMD))
  
  return(metrics)
}