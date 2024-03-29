#' Function to Identify Network Regulators from Undirected Weighted Networks
#' 
#' Identifies network un-weighted regulators from Weighted undirected networks
#' 
#' @param adj Required. An n x n weighted upper triangular adjacency in the matrix 
#' class format.
#' @param G Required. A named vector of node scores.
#' @param h Optional. Neighborhood search distance (h nodes away from current node) 
#' (Default = 3)
#' @param n Optional. number of permuations to run. (Default = 100)
#' @param FDR Optional. Adjusted pvalue cutoff for regulator selection.
#' (Default = 0.05)
#' 
#' @return scores = n x 4 dimensional list with columns giving neighborhood
#'based score, adjusted pvalue, whether a gene is regulator/global regulator.
#'
#' @importFrom foreach %dopar%
#' @export
regulatorAnalysis.undirected_weighted <- function(adj, G, h = 3, n = 100, FDR = 0.05){
  
  # Convert adjacency to igraph object
  g = igraph::graph_from_adjacency_matrix(adj, mode = 'upper')
  background.genes = igraph::V(g)$name
  
  # Compute node scores based on neighborhood scores
  node.scores = score.nodes(g, G, h)
  
  # Permute node labels and calculate null scores
  i <- NULL
  perm.node.scores = foreach::foreach(i = 1:n, .combine = cbind, .packages = c('igraph', 'dplyr', 'parallel', 'doParallel', 'foreach')) %dopar% {
    pg = igraph::permute(g, sample(1:length(background.genes), length(background.genes)))
    perm.node.scores = score.nodes(pg, G, h)
  } 
  i <- NULL
  # Perform one sample t-test to estimate significance
  pval = foreach::foreach(i = 1:n, .combine = rbind) %dopar% {
    tmp = stats::t.test(perm.node.scores[i,], mu = node.scores[i], alternative = 'less')
    data.frame(pval = tmp$p.value, t = tmp$statistic, t.low = tmp$conf.int[1], t.high = tmp$conf.int[2])
  }
  fdr = stats::p.adjust(pval$pval, method = 'fdr')
  names(fdr) = background.genes
  
  # Calculate node degree for identifying global regulators
  node.degree = igraph::degree(g)
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = stats::sd(node.degree, na.rm = T)
  
  # Promote high degree nodes as global regulators 
  key.regulators = list()
  key.regulators$scores = node.scores
  key.regulators$fdr = fdr
  key.regulators$global.regulators = background.genes[(node.degree > (mean.node.degree + 2*stddev.node.degree)) & (fdr <= FDR)]
  key.regulators$local.regulators = background.genes[(node.degree < (mean.node.degree + 2*stddev.node.degree)) & (fdr <= FDR)]
  
  return(key.regulators)
}