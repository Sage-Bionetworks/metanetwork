#' Function to Identify Network Regulators from Directed Networks
#' 
#' Identifies network un-weighted regulators from directed networks
#' 
#' @param adj Required. An n x n weighted upper triangular adjacency in the matrix 
#' class format.
#' @param G Required. A named vector of node scores.
#' @param h Optional. Neighborhood search distance (h nodes away from current node) 
#' (Default = 3)
#' @param FDR Optional. Adjusted pvalue cutoff for regulator selection.
#' (Default = 0.05)
#' 
#' @return scores = n x 5 dimensional data frame with columns giving neighborhood
#'based score, adjusted pvalue, whether a gene is regulator/global regulator.
#'
#' @importFrom foreach %dopar%
#' @export
regulatorAnalysis.directed <- function(adj, G, h=3, FDR = 0.05){
  
  # Convert adjacency to igraph object
  g = igraph::graph_from_adjacency_matrix(adj, mode = 'directed')
  background.genes = igraph::V(g)$name
  
  # Find H-Layer neighborhood
  neighbor.nodes = lapply(1:h, function(hi,sg){
    igraph::ego(sg, order = hi, nodes = igraph::V(sg), mode = 'out')
  }, g)
  
  # Perform enrichment analysis for every gene in every layer and pick the minimum p-value
  fdr = sapply(1:length(neighbor.nodes), function(x, neighborNodes, G, backgroundGenes){
    foreach::foreach(i = 1:length(neighborNodes[[x]]), 
                     .combine = c, 
                     .export =c('neighborNodes', 'G', 'backgroundGenes', 'fisherEnrichment')) %dopar% {
                       metanetwork::fisherEnrichment(names(neighborNodes[[x]][i][[1]]), G, backgroundGenes)$pval
                     }
  }, neighbor.nodes, G, background.genes) %>%
    apply(1, min, na.rm = T) %>%
    stats::p.adjust(method = 'fdr')
  names(fdr) = background.genes
  
  # Calculate node degree for identifying global regulators
  node.degree = igraph::degree(g)
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = stats::sd(node.degree, na.rm = T)
  
  # Find leaf nodes to identify global regulators
  node.in.degree = igraph::degree(g, mode = 'in')
  
  key.regulators = list()
  key.regulators$fdr = fdr
  key.regulators$regulators = background.genes[(fdr <= FDR)]
  key.regulators$global.regulators = background.genes[((node.degree > (mean.node.degree + 2*stddev.node.degree)) | (node.in.degree == 0)) & (fdr <= FDR)]
  key.regulators$local.regulators = setdiff(key.regulators$regulators, key.regulators$global.regulators)
  
  return(key.regulators)
}