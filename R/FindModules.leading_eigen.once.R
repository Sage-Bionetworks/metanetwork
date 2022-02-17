#' Find Modules with Leading Eigen Edges
#' 
#' This function tries to find densely connected subgraphs in a graph by 
#' calculating the leading non-negative eigenvector of the modularity matrix of 
#' the graph.
#' 
#' @inheritParams findModules.edge_betweenness.once
#' 
#' @return  GeneModules = n x 3 dimensional data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
#' 
#' @importFrom magrittr %>%
#' @export
findModules.leading_eigen.once <- function(adj, min.module.size){
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'undirected', weighted = T, diag = F)
  
  # Get connected graph (remove unconnected nodes)
  node.degree = igraph::degree(g)
  unconnected.genes = data.frame(Gene.ID = igraph::V(g)$name[node.degree == 0],
                                 moduleNumber = 0)
  
  # Get modules using leading_eigen algorithm (MEJ Newman: Finding community structure using the eigenvectors of matrices, Physical Review E 74 036104, 2006)
  sg = induced_subgraph(g, which(node.degree != 0))
  mod = igraph::cluster_leading_eigen(sg)
  
  # Get individual clusters from the igraph community object
  geneModules = igraph::membership(mod) %>%
    unclass %>%
    as.data.frame %>%
    plyr::rename(c('.' = 'moduleNumber'))
  
  geneModules = cbind(data.frame(Gene.ID = rownames(geneModules)),
                      geneModules)              
  geneModules = rbind(geneModules, unconnected.genes)
  
  # Rename modules with size less than min module size to 0
  filteredModules = geneModules %>% 
    dplyr::group_by(moduleNumber) %>%
    dplyr::summarise(counts = length(unique(Gene.ID))) %>%
    dplyr::filter(counts >= min.module.size)
  geneModules$moduleNumber[!(geneModules$moduleNumber %in% filteredModules$moduleNumber)] = 0
  
  # Change cluster number to color labels
  geneModules$moduleLabel = WGCNA::labels2colors(geneModules$moduleNumber)
  
  return(geneModules)
}