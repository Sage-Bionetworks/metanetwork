#' Find Modules with Network Adjacency Matrix Using spinglass Clustering
#' 
#' This function tries to get modules from network adjacency matrix using igraph's 
#' spinglass clusting function.
#' 
#' @inheritParams findModules.edge_betweenness.once
#' 
#' @return  GeneModules = n x 3 dimensional data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
#' 
#' @importFrom magrittr %>%
#' @export
findModules.spinglass.once <- function(adj, min.module.size){
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'undirected', weighted = T, diag = F)
  
  # Find connected components
  scc = igraph::components(g)
  
  # Find modules for each component using spinglass algorithm (http://arxiv.org/abs/cond-mat/0603718)
  mod = lapply(unique(scc$membership), function(x, g, scc){
    sg = igraph::induced_subgraph(g, which(scc$membership == x))
    if (sum(scc$membership == x) == 1){
      geneModules = data.frame(Gene.ID = igraph::V(g)$name[scc$membership == x],
                               moduleNumber = 0)
    } else{
      mod = igraph::cluster_spinglass(sg)
      geneModules = data.frame(Gene.ID = igraph::V(sg)$name,
                               moduleNumber = unclass(igraph::membership(mod)))
    }
  }, g, scc)
  
  mod.sz = c(0, cumsum(sapply(mod, function(x) max(x$moduleNumber))))
  mod.sz = mod.sz[1:(length(mod.sz) -1)]
  geneModules = mapply(function(x,y){
    x$moduleNumber = x$moduleNumber + y
    return(x)
  }, mod, mod.sz, SIMPLIFY = F) %>%
    rbindlist(use.names = T, fill = T)
  
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