key.driver.analysis <- function(adj, G, h=3, FDR = 0.05){
  
  key.drivers = list()
  # Convert adjacency to igraph object
  g = igraph::graph_from_adjacency_matrix(adj, mode = 'upper')
  background.genes = V(g)$name
  
  # Find H-Layer neighborhood
  neighbor.nodes = lapply(1:h, function(hi,sg){
    ego(sg, order = hi, nodes = V(sg), mode = 'all')
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
  
  # Calculate node degree for identifying global drivers
  node.degree = degree(g)
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = sd(node.degree, na.rm = T)
  
  key.drivers$fdr = p.val
  key.drivers$global.regulators = names(p.val)[(node.degree > (mean.node.degree + 2*stddev.node.degree)) & (p.val <= FDR)]
  key.drivers$local.regulators = names(p.val)[(node.degree < (mean.node.degree + 2*stddev.node.degree)) & (p.val <= FDR)]
  
  return(key.drivers)
}