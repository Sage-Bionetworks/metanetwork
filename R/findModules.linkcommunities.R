# Function to get modules from network adjacency matrix using Link communities algorithm
findModules.linkcommunities <- function(adj, min.module.size = 20){
  
  # Note: For this function to work get the package linkcomm from CRAN
  
  # Input
  #      adj = n x n adjacency matrix in the ltCMatrix format
  #      min.module.size = integer specifying the minimum module size
  
  # Output (list of following elements)
  #      modules = data frame of dimension n x 3 with columns named GeneIDs, moduleNumber and moduleLabel
  
  # Error functions
  if(class(adj) != "ltCMatrix")
    stop('Adjacency matrix should be of class ltCMatrix')
  
  if(dim(adj)[1] != dim(adj)[2])
    stop('Adjacency matrix should be symmetric')
  
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'upper', weighted = NULL, diag = F)
  
  ## Get edgelist
  elist = igraph::as_edgelist(g)
  
  ## Run link communities function
  comm = linkcomm::getLinkCommunities(elist)
  
  # Get individual clusters from the community object
  geneModules = comm$nodeclusters %>%
    dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
    plyr::rename(c('node' = 'Gene.ID', 'cluster' = 'moduleNumber'))

  # Add missing genes
  Gene.ID = setdiff(igraph::V(g)$name, geneModules$Gene.ID)
  geneModules = rbind(geneModules, 
                      data.frame(Gene.ID = Gene.ID,
                                 moduleNumber = max(geneModules$moduleNumber, na.rm = T) + seq(1,length(Gene.ID))))
  
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