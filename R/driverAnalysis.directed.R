driverAnalysis.directed <- function(g, G, h=3, FDR = 0.05){
  
  # Convert adjacency to igraph object
  background.genes = igraph::V(g)$name
  
  # Find H-Layer neighborhood
  neighbor.nodes = lapply(1:h, function(hi,sg){
    igraph::ego(sg, order = hi, nodes = igraph::V(sg), mode = 'out')
  },g)
  
  # Perform enrichment analysis for every gene in every layer and pick the minimum p-value
  fdr = sapply(1:length(neighbor.nodes), function(x, neighborNodes, G, backgroundGenes){
    foreach::foreach(i = 1:length(neighborNodes[[x]]), .combine = c, .export =c('neighborNodes', 'G', 'backgroundGenes', 'fisherEnrichment')) %dopar% {
      fisherEnrichment(names(neighborNodes[[x]][i][[1]]), G, backgroundGenes)$pval
    }
  }, neighbor.nodes, G, background.genes) %>%
    apply(1, min, na.rm = T) %>%
    stats::p.adjust(method = 'fdr')
  names(fdr) = background.genes
  
  # Calculate node degree for identifying global drivers
  node.degree = igraph::degree(g)
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = sd(node.degree, na.rm = T)
  
  # Find leaf nodes to identify global drivers
  node.in.degree = degree(g, mode = 'in')
  
  key.drivers = list()
  key.drivers$fdr = fdr
  key.drivers$regulators = background.genes[(fdr <= FDR)]
  key.drivers$global.regulators = background.genes[((node.degree > (mean.node.degree + 2*stddev.node.degree)) | (node.in.degree == 0)) & (fdr <= FDR)]
  key.drivers$local.regulators = setdiff(key.drivers$regulators, key.drivers$global.regulators)
  
  return(key.drivers)
}