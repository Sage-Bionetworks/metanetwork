# Function to calculate modules
getModules <- function(adjMat, 
                       TOMType = 'unsigned', #other options are 'signed'
                       TOMDenom = 'min', #other options are 'mean',
                       linkageType = 'ward', #
                       distanceType = 'euclidean', #
                       minClusterSize = 30,
                       deepSplit = F){
  
  # Check adjMat is of class matrix
  if (!is(adjMat, 'matrix'))
    stop('adjMat should be of class matrix')
  
  # Get topological overlap matrix
  TOM = WGCNA::TOMdist(adjMat, 
                       TOMType, 
                       TOMDenom, 
                       verbose = 2)
  colnames(TOM) = colnames(adjMat)
  rownames(TOM) = rownames(adjMat)
  
  # Get hierarchical clusters from TOM
  collectGarbage()
  Rclusterpp::Rclusterpp.setThreads(threads=4)
  clust = Rclusterpp::Rclusterpp.hclust(TOM,
                                        method = linkageType,
                                        distance = distanceType)
  
  # Get individual clusters from the hierarchical tree
  collectGarbage()
  clust.numLabels = try(dynamicTreeCut::cutreeDynamic(clust,
                                                      minClusterSize = minClusterSize,
                                                      method = 'tree',
                                                      deepSplit = deepSplit))
  
  if (is(clust.numLabels,'try-error'))
    clust.numLabels = rep(0, dim(adjMat)[1])
  
  # Change cluster number to color labels
  collectGarbage()
  labels = WGCNA::labels2colors(clust.numLabels)
  
  # Get results
  geneModules = data.frame(GeneIDs = rownames(adjMat),
                           moduleNumber = clust.numLabels, 
                           modulelabels = labels,
                           row.names = rownames(adjMat))
  
  return(list(geneModules = geneModules,
              TOM = TOM,
              clust = clust))
}