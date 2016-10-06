regulatorAnalysis.directed_weighted <- function(g, G, h = 3, n = 100, FDR = 0.05){
  
  # Convert adjacency to igraph object
  background.genes = igraph::V(g)$name
  
  # Compute node scores based on neighborhood scores
  node.scores = score.nodes(g, G, h, mode = 'out')
  
  # Permute node labels and calculate null scores
  perm.node.scores = foreach::foreach(i = 1:n, 
                                      .combine = cbind, 
                                      .packages = c('igraph', 'dplyr', 'parallel', 'doParallel', 'foreach'),
                                      .verbose = TRUE) %dopar% {
                                        pg = igraph::permute(g, sample(1:length(background.genes), length(background.genes)))
                                        perm.node.scores = score.nodes(pg, G, h, mode = 'out')
                                      } 
  
  # Perform one sample t-test to estimate significance
  pval = foreach::foreach(i = 1:length(background.genes), .combine = rbind) %dopar% {
    tmp = t.test(perm.node.scores[i,], mu = node.scores[i], alternative = 'less')
    data.frame(pval = tmp$p.value, t = tmp$statistic, t.low = tmp$conf.int[1], t.high = tmp$conf.int[2])
  }
  fdr = p.adjust(pval$pval, method = 'fdr')
  names(fdr) = background.genes
  
  # Calculate node degree for identifying global regulators
  node.degree = igraph::degree(g)
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = sd(node.degree, na.rm = T)
  
  # Find leaf nodes to identify global regulators
  node.in.degree = degree(g, mode = 'in')
  
  # Promote high degree nodes as global regulators 
  key.regulators = list()
  key.regulators$scores = node.scores
  key.regulators$fdr = fdr
  key.regulators$regulators = background.genes[(fdr <= FDR)]
  key.regulators$global.regulators = background.genes[((node.degree > (mean.node.degree + 2*stddev.node.degree)) | (node.in.degree == 0)) & (fdr <= FDR)]
  key.regulators$local.regulators = setdiff(key.regulators$regulators, key.regulators$global.regulators)
  
  return(key.regulators)
}