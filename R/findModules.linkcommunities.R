# Function to get modules from network adjacency matrix using Link communities algorithm
findModules.linkcommunities <- function(adj, min.module.size = 3){
  
  # Note: For this function to work get the package linkcomm from CRAN
  
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
  
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'upper', weighted = T, diag = F)
  
  ## Get edgelist
  elist = igraph::as_edgelist(g)
  
  ## Run link communities function
  comm = linkcomm::getLinkCommunities(elist)
  
  ## Get edge clusters
  eclust = cutree(comm$hclust, h = comm$pdmax)
  names(eclust) = igraph::E(comm$igraph)
  
  ## Get node communities
  nodes = lapply(seq(max(eclust)), function(x, eclust, comm){
    tmp.g = igraph::subgraph.edges(comm$igraph, eids = which(eclust == x), delete.vertices = T)
    data.frame(Gene.ID = igraph::V(tmp.g)$name) %>%
      dplyr::mutate(moduleNumber = x,
                    moduleSize = length(unique(Gene.ID)))
  }, eclust, comm) %>%
    data.table::rbindlist(use.names = T, fill = T)
  
  ## Combine smaller modules to form a no module
  nodes$moduleNumber[nodes$moduleSize < min.module.size] = 0
  
  # Get individual clusters from the community object
  geneModules = dplyr::filter(nodes, moduleNumber != 0) %>%
    dplyr::group_by(Gene.ID) %>%
    dplyr::top_n(1, moduleSize) %>%
    dplyr::top_n(1, moduleNumber) %>%
    dplyr::select(-moduleSize) %>%
    dplyr::mutate(moduleNumber = factor(moduleNumber),
                  moduleNumber = as.numeric(moduleNumber))
  
  # Add missing genes
  Gene.ID = setdiff(igraph::V(g)$name, geneModules$Gene.ID)
  geneModules = rbind(data.frame(geneModules), 
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