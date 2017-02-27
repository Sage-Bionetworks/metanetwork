# Function to get modules from network adjacency matrix
findModules.hclust <- function(adj, agglom.method = 'ward', clustDistance = 'euclidean', min.module.size = 3){
  # Input
  #      adj = n x n upper triangular adjacency in the matrix class format
  #      min.module.size = integer between 1 and n genes 
  
  # Output
  #      geneModules = n x 3 dimensional data frame with column names as Gene.ID, moduleNumber, and moduleLabel
  
  # Error functions
  if(class(adj) != "matrix")
    stop('Adjacency matrix should be of class matrix')
  
  if(dim(adj)[1] != dim(adj)[2])
    stop('Adjacency matrix should be symmetric')
  
  if(!all(adj[lower.tri(adj)] == 0))
    stop('Adjacency matrix should be upper triangular')
  
  adj = adj + t(adj)
  
  TOM = WGCNA::TOMsimilarity(adj);
  dissTOM = 1-TOM
  
  dissStruct = dist(dissTOM, method = clustDistance)
  
  geneTree = flashClust::hclust(dissStruct, method = agglom.method)
  
  mod = dynamicTreeCut::cutreeDynamic(dendro = geneTree, 
                                          method = 'hybrid',
                                          distM = dissTOM,
                                          pamRespectsDendro = F,
                                          minClusterSize = min.module.size)
  names(mod) = rownames(adj)
  geneModules = data.frame(Gene.ID = names(mod),
                           moduleNumber = mod)
  
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