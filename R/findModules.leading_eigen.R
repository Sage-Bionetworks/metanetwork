# Function to get modules from network adjacency matrix
findModules.leading_eigen <- function(adj, min.module.size = 20){
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'upper', weighted = NULL, diag = F)
  
  # Get connected graph (remove unconnected nodes)
  node.degree = igraph::degree(g)
  unconnected.genes = data.frame(Gene.ID = igraph::V(g)$name[node.degree == 0],
                                 moduleNumber = 0)
  
  # Get modules using leading_eigen algorithm (MEJ Newman: Finding community structure using the eigenvectors of matrices, Physical Review E 74 036104, 2006)
  sg = induced_subgraph(g, which(node.degree != 0))
  mod = igraph::cluster_leading_eigen(sg)
  
  # Get individual clusters from the igraph community object
  geneModules = igraph::membership(mod) %>%
    unclass %>%
    as.data.frame %>%
    plyr::rename(c('.' = 'moduleNumber'))
  
  geneModules = cbind(data.frame(Gene.ID = rownames(geneModules)),
                      geneModules)              
  geneModules = rbind(geneModules, unconnected.genes)
  
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