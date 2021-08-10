# Function to get modules from network adjacency matrix
findModules.megena <- function(adj, method = "pearson", FDR.cutoff = 0.05, module.pval = 0.05, hub.pval = 0.05,doPar = TRUE){
  # Input
  #      adj = n x n upper triangular adjacency in the matrix class format
  #      method = method for correlation. either pearson or spearman. 
  #      FDR.cutoff = FDR threshold to define significant correlations upon shuffling samples.
  #      module.pval = module significance p-value. Recommended is 0.05. 
  #      hub.pval = connectivity significance p-value based random tetrahedral networks
  #      doPar = if to allow parallelization of clusters
  # Output
  #      geneModules = n x 3 dimensional data frame with column names as Gene.ID, moduleNumber, and moduleLabel
  library(MEGENA)
  library(igraph)
  cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
  hub.perm = 100; # number of permutations for calculating connectivity significance p-value. 
  
  g = igraph::graph.adjacency(adj, mode = 'undirected', weighted = T, diag = F)
  
  MEGENA.output <- do.MEGENA(g,
                             mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = FALSE,
                             min.size = 30,max.size = vcount(g)/2,
                             doPar = TRUE,num.cores = nc,n.perm = hub.perm,
                             save.output = FALSE)

  geneModules <- module_convert_to_table(MEGENA.output,mod.pval = 0.05,
                          hub.pval = 0.05,min.size = 30,max.size=vcount(g)/2)
  count_mtx <- table(geneModules['module.parent'])
  colnames(count_mtx) <- c('moduleNumber','moduleSize')
  colnames(geneModules) <- c('Gene.ID','moduleNumber','nodeDegree','nodeStrength','hub','module')
  geneModules['moduleSize'] <- 0
  for (i in 1:nrow(geneModules)){
    temp = as.character(geneModules$moduleNumber[i])
    geneModules$moduleSize = as.integer(count_mtx[temp])
  }
  
  geneModules = geneModules %>%
    group_by(Gene.ID) %>%
    dplyr::top_n(1, moduleSize) %>%
    dplyr::top_n(1, moduleNumber) %>%
    dplyr::select(-moduleSize) %>%
    dplyr::mutate(moduleNumber = factor(moduleNumber),
                  moduleNumber = as.numeric(moduleNumber))
  
  
  # Rename modules with size less than min module size to 0
  filteredModules = geneModules %>% 
    dplyr::group_by(moduleNumber) %>%
    dplyr::summarise(counts = length(unique(Gene.ID))) %>%
    dplyr::filter(counts >= 30)
  geneModules$moduleNumber[!(geneModules$moduleNumber %in% filteredModules$moduleNumber)] = 0
  
  # Change cluster number to color labels
  geneModules$moduleNumber = as.numeric(factor(geneModules$moduleNumber))
  geneModules$moduleLabel = WGCNA::labels2colors(geneModules$moduleNumber)
  mod = geneModules[c('Gene.ID','moduleNumber','moduleLabel')]

  return(mod)
}