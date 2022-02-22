#' Find Modules with Network Adjacency Single iteration
#' 
#' This function finds modules from network adjacency matrix.
#' 
#' @param adj A n x n upper triangular adjacency in the matrix class format.
#' @param min.module.size Optional. Integer between 1 and n genes. (Default = 30)
#'
#' @return  GeneModules = n x 3 dimensional data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
#' 
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
findModules.edge_betweenness.once <- function(adj, min.module.size){
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'undirected', weighted = T, diag = F)
  
  # Get modules using edge betweenness algorithm (M Newman and M Girvan: Finding and evaluating community structure in networks, Physical Review E 69, 026113 (2004))
  mod = igraph::cluster_edge_betweenness(g)
  
  # Get individual clusters from the igraph community object
  geneModules = igraph::membership(mod) %>%
    unclass %>%
    as.data.frame %>%
    plyr::rename(c('.' = 'moduleNumber'))
  
  geneModules = cbind(data.frame(Gene.ID = rownames(geneModules)),
                      geneModules)              
  
  # Rename modules with size less than min module size to 0
  filteredModules = geneModules %>% 
    dplyr::group_by(.data$moduleNumber) %>%
    dplyr::summarise(counts = length(unique(.data$Gene.ID))) %>%
    dplyr::filter(.data$counts >= min.module.size)
  geneModules$moduleNumber[!(geneModules$moduleNumber %in% filteredModules$moduleNumber)] = 0
  
  # Change cluster number to color labels
  geneModules$moduleLabel = WGCNA::labels2colors(geneModules$moduleNumber)
  
  return(geneModules)
}