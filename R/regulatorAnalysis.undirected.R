#' Function to Identify Network Regulators from Undirected Networks
#' 
#' Identifies network un-weighted regulators from Undirected networks
#' 
#' @param adj Required. An n x n weighted upper triangular adjacency in the matrix 
#' class format.
#' @param G Required. A named vector of node scores.
#' @param h Optional. Neighborhood search distance (h nodes away from current node) 
#' (Default = 3)
#' @param FDR Optional. Adjusted pvalue cutoff for regulator selection.
#' (Default = 0.05)
#' 
#' @return scores = n x 3 dimensional list with columns giving neighborhood
#'based score, adjusted pvalue, whether a gene is regulator/global regulator.
#'
#' @importFrom foreach %dopar%
#' @export
regulatorAnalysis.undirected <- function(adj, G, h=3, FDR = 0.05){
  
  # Convert adjacency to igraph object
  g = igraph::graph_from_adjacency_matrix(adj, mode = 'upper')
  background.genes = igraph::V(g)$name
  
  # Find H-Layer neighborhood
  neighbor.nodes = lapply(1:h, function(hi,sg){
    igraph::ego(sg, order = hi, nodes = igraph::V(sg), mode = 'all')
  },g)
  
  # Perform enrichment analysis for every gene in every layer and pick the minimum p-value
  p.val = foreach (x = 1:length(neighbor.nodes)) %dopar% {
    p.val = sapply(neighbor.nodes[[x]], function(y, G, backgroundGenes){
      metanetwork::fisherEnrichment(names(y), G, backgroundGenes)$pval
    }, G, background.genes)
  } %>%
    do.call(cbind,.) %>%
    apply(1, min, na.rm = T) %>%
    stats::p.adjust(method = 'fdr')
  names(p.val) = background.genes
  
  # Calculate node degree for identifying global regulators
  node.degree = igraph::degree(g)
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = sd(node.degree, na.rm = T)
  
  # Coallate results
  key.regulators = list()
  key.regulators$fdr = p.val
  key.regulators$global.regulators = names(p.val)[(node.degree > (mean.node.degree + 2*stddev.node.degree)) & (p.val <= FDR)]
  key.regulators$local.regulators = names(p.val)[(node.degree < (mean.node.degree + 2*stddev.node.degree)) & (p.val <= FDR)]
  
  return(key.regulators)
}