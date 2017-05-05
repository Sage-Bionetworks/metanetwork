# Function to identify network regulators from directed weighted networks
regulatorAnalysis.directed_weighted <- function(g, G, h = 3, n = 100, correction.method = 'bonferroni', pval.cutoff = 0.01){
  
  # Input
  #      adj = n x n weighted upper triangular adjacency in the matrix class format
  #      G = a named vector of node scores 
  #      h = neighborhood search distance (h nodes away from current node)
  #      n = number of randomisation for pvalue computation
  #      correction.method = multiple testing correction method
  #      pval.cutoff = adjusted pvalue cutoff for regulator selection
  
  # Output
  #      scores = n x 5 dimensional data frame with columns giving neighborhood based score, adjusted pvalue, whether a gene is regulator/global regulator
  
  # Error checking
  if(dim(adj)[1] != dim(adj)[2])
    stop('Adjacency matrix should be symmetric')
  
  if(!all(adj[lower.tri(adj)] == 0))
    stop('Adjacency matrix should be upper triangular')
  
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'upper', weighted = T, diag = F)
  
  # Get the background genes
  background.genes = igraph::V(g)$name
  
  # Compute node scores based on neighborhood scores
  node.scores = score.nodes(g, G, h, mode = 'out')
  
  # Permute node labels and calculate null distribution of scores
  perm.node.scores = foreach::foreach(i = 1:n,
                                      .combine = cbind, 
                                      .packages = c('igraph', 'dplyr', 'parallel', 'doParallel', 'foreach'),
                                      .verbose = FALSE) %dopar% {
                                        pg = igraph::permute(g, sample(1:length(background.genes), length(background.genes)))
                                        perm.node.scores = score.nodes(pg, G, h, mode = 'out')
                                      } 
  
  # Perform one sample t-test to estimate significance
  pval = foreach::foreach(i = 1:length(background.genes), .combine = rbind) %dopar% {
    tmp = t.test(perm.node.scores[i,], mu = node.scores[i], alternative = 'less')
    data.frame(pval = tmp$p.value, t = tmp$statistic, t.low = tmp$conf.int[1], t.high = tmp$conf.int[2])
  }
  adj.P.Val = p.adjust(pval$pval, method = correction.method)
  names(adj.P.Val) = background.genes
  
  scores = data.frame(Gene.ID = names(node.scores),
                      scores = node.scores,
                      adj.P.Val = adj.P.Val[names(node.scores)])
  scores$regulator = FALSE
  scores$regulator[scores$adj.P.Val <= pval.cutoff] = TRUE
  
  # Calculate node degree for identifying global regulators
  node.degree = igraph::strength(g, mode = 'all')
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = sd(node.degree, na.rm = T)
  
  # Find leaf nodes to identify global regulators
  node.in.degree = igraph::degree(g, mode = 'in')
  
  # Promote high degree nodes as global regulators 
  scores$global.regulator = FALSE
  scores$global.regulator[scores$regulator == TRUE & (node.in.degree == 0)] = TRUE
  scores$global.regulator[scores$regulator == TRUE & 
                             (node.degree > (mean.node.degree + 2*stddev.node.degree))] = TRUE
  
  return(dplyr::arrange(scores, desc(global.regulator), desc(regulator), adj.P.Val))
}