#' Find Modules with Network Adjacency Matrix Using Link Communities Algorithm
#' 
#' This function tries to get modules from network adjacency matrix using Link 
#' communities algorithm.
#' 
#' @inheritParams findModules.edge_betweenness.once
#' 
#' @return  GeneModules = n x 3 dimensional data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
#' 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
findModules.linkcommunities.once <- function(adj, min.module.size){
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'undirected', weighted = T, diag = F)
  
  ## Get edgelist
  elist = igraph::as_edgelist(g)
  
  ## Run link communities function
  comm = linkcomm::getLinkCommunities(elist)
  
  ## Get edge clusters
  eclust = stats::cutree(comm$hclust, h = comm$pdmax)
  names(eclust) = igraph::E(comm$igraph)
  
  ## Filter communities less than min.module.size
  comm.to.remove = which(table(eclust) <= min.module.size)
  eclust[eclust %in% comm.to.remove] = 0
  
  nodes = lapply(setdiff(unique(eclust), 0), function(x, eclust, comm){
    tmp.g = igraph::subgraph.edges(comm$igraph, eids = which(eclust == x), delete.vertices = T)
    data.frame(Gene.ID = igraph::V(tmp.g)$name) %>%
      dplyr::mutate(moduleNumber = x,
                    moduleSize = length(unique(.data$Gene.ID)))
  }, eclust, comm) %>%
    data.table::rbindlist(use.names = T, fill = T)
  
  ## Combine smaller modules to form a no module
  nodes$moduleNumber[nodes$moduleSize < min.module.size] = 0
  
  # Get individual clusters from the community object
  geneModules = dplyr::filter(nodes, .data$moduleNumber != 0) %>%
    dplyr::group_by(Gene.ID) %>%
    dplyr::top_n(1, .data$moduleSize) %>%
    dplyr::top_n(1, .data$moduleNumber) %>%
    dplyr::select(-.data$moduleSize) %>%
    dplyr::mutate(moduleNumber = factor(moduleNumber),
                  moduleNumber = as.numeric(moduleNumber))
  
  # Add missing genes
  Gene.ID = setdiff(igraph::V(g)$name, geneModules$Gene.ID)
  geneModules = rbind(data.frame(geneModules), 
                      data.frame(Gene.ID = Gene.ID,
                                 moduleNumber = 0))
  
  # Rename modules with size less than min module size to 0
  filteredModules = geneModules %>% 
    dplyr::group_by(.data$moduleNumber) %>%
    dplyr::summarise(counts = length(unique(.data$Gene.ID))) %>%
    dplyr::filter(.data$counts >= min.module.size)
  geneModules$moduleNumber[!(geneModules$moduleNumber %in% filteredModules$moduleNumber)] = 0
  
  # Change cluster number to color labels
  geneModules$moduleNumber = as.numeric(factor(geneModules$moduleNumber))
  geneModules$moduleLabel = WGCNA::labels2colors(geneModules$moduleNumber)
  
  return(geneModules)
}