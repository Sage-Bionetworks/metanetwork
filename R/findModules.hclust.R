#' Find Modules with Hierarchical Clustering
#' 
#' This function finds modules from Hierarchical Clustering `fastcluster::hclust.vector()`
#' 
#' @param adj Required. A n x n upper triangular adjacency in the matrix class format.
#' @param aggloMethod Optional. The agglomeration method to be used. This must be 
#' (an unambiguous abbreviation of) one of "single", "ward", "centroid" or "median".
#' (Default = "ward")
#' @param clustDistance Optional. The distance measure to be used. This must be 
#' one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". 
#' Any unambiguous substring can be given. (Default = "euclidean")
#' @param minModuleSize Optional. Integer between 1 and n genes. (Default = 30)
#'
#' @return  GeneModules = n x 3 dimensional data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
#' 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
findModules.hclust <- function(adj, aggloMethod = 'ward', clustDistance = 'euclidean', minModuleSize = 30){
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
  
  # Use fast hierarchichal clustering
  geneTree = fastcluster::hclust.vector(adj, method = aggloMethod, metric = clustDistance)
  
  # Find distance between matrix
  distAdj = stats::dist(adj, method = clustDistance, diag = T, upper = T)
  
  # Cut tree to form clusters
  mod = dynamicTreeCut::cutreeDynamic(dendro = geneTree, 
                                      minClusterSize = minModuleSize,
                                      method = 'hybrid',
                                      distM = as.matrix(distAdj),
                                      deepSplit =FALSE,
                                      pamRespectsDendro = FALSE)
  
  geneModules = data.frame(Gene.ID = rownames(adj),
                           moduleNumber = mod)
  
  # Rename modules with size less than min module size to 0
  filteredModules = geneModules %>% 
    dplyr::group_by(.data$moduleNumber) %>%
    dplyr::summarise(counts = length(unique(.data$Gene.ID))) %>%
    dplyr::filter(.data$counts >= minModuleSize)
  geneModules$moduleNumber[!(geneModules$moduleNumber %in% filteredModules$moduleNumber)] = 0

  # Change cluster number to color labels
  geneModules$moduleLabel = WGCNA::labels2colors(geneModules$moduleNumber)
  
  return(geneModules)
}