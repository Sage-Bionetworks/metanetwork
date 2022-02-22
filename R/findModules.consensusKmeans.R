#' Finds Moduless With Kmeans Clustering
#' 
#' Function to get consensus modules from individual partition matrices
#' 
#' @param partition.adj Required. A partition.adj = n x m adjacency matrix, where 
#' n is the number of genes and m = number of clustering methods * number of 
#' clusters in each method.
#' @param min.module.size Optional. An integer between 1 and n genes (Default = 20)
#' @param usepam Optional. A logical for input into pam based kmeans clustering 
#' to find the number of clusters with the function `fpc::pamk`. If TRUE, pam is
#' used, otherwise clara (recommended for large datasets with 2,000 or more 
#' observations; dissimilarity matrices can not be used with clara). (Default = 20)
#' 
#' @return A dataframe of Gene Modules
#' 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
findModules.consensusKmeans <- function(partition.adj, min.module.size = 20, usepam = FALSE){
  # Input
  #      partition.adj = n x m adjacency matrix, where n is the number of genes and m = number of clustering methods * number of clusters in each method
  #      min.module.size = integer between 1 and n genes 
  
  # Output
  #      geneModules = n x 3 dimensional data frame with column names as Gene.ID, moduleNumber, and moduleLabel
  
  # Error functions
  if(class(partition.adj) != "matrix")
    stop('partition.adjacency matrix should be of class matrix')
  
  # Use pam based kmeans clustering to find the number of clusters
  mod = fpc::pamk(partition.adj, krange = 2:30, usepam = usepam)
  
  # Get individual clusters from the igraph community object
  geneModules = data.frame(Gene.ID = names(mod$pamobject$cluster),
                           moduleNumber = mod$pamobject$cluster)              
  
  # Rename modules with size less than min module size to 0
  filteredModules = geneModules %>% 
    dplyr::group_by(.data$moduleNumber) %>%
    dplyr::summarise(counts = length(unique(.data$Gene.ID))) %>%
    dplyr::filter(.data$counts >= min.module.size)
  geneModules$moduleNumber[!(geneModules$moduleNumber %in% filteredModules$moduleNumber)] = 0
  
  # Change cluster number to color labels
  geneModules$moduleLabel = WGCNA::labels2colors(geneModules$moduleNumber)
  
  return(geneModules)
}