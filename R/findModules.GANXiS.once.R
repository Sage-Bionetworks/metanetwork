#' Find Modules using GANXIS
#' 
#' This function finds modules with GANXIS
#' 
#' @inheritParams findModules.edge_betweenness.once
#' 
#' @return  GeneModules = n x 3 dimensional data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
#' 
#' @importFrom magrittr %>%
#' @export
findModules.GANXiS.once <- function(adj, path, min.module.size){
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'undirected', weighted = T, diag = F)
  
  # Get modules using GANXiS
  system('rm -rf ./tmp')
  system('mkdir ./tmp')
  
  ## Write network as an edgelist to a file
  elist = igraph::as_edgelist(g)
  elist = cbind(elist, 1)
  utils::write.table(elist, file = './tmp/input.txt', row.names = F, col.names = F, quote=F)
  
  ## Run GANXiS
  system(paste('java','-jar',paste0(path, 'GANXiSw.jar'),
               '-i','./tmp/input.txt',
               '-d', './tmp',
               '-Sym','1',
               '-ov','0',
               '-Onc','1',
               '-seed', '123456789'))
  
  ## Get output 
  mod = utils::read.table('./tmp/SLPAw_input_run1_r0.5_v3_T100.icpm.node-com.txt')
  system('rm -rf ./tmp')
  
  # Get individual clusters from the igraph community object
  numLabels = mod$V2+1
  names(numLabels) = mod$V1
  
  # Get individual clusters from the igraph community object
  geneModules = numLabels %>%
    unclass %>%
    as.data.frame %>%
    plyr::rename(c('.' = 'moduleNumber'))
  
  geneModules = cbind(data.frame(Gene.ID = rownames(geneModules)),
                      geneModules)              
  
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