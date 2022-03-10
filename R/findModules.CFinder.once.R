#' This function tries to get modules from network adjacency matrix using Link 
#' communities algorithm.
#' 
#' @param i Required. Character name of temp file name to generate.
#' @inheritParams findModules.edge_betweenness.once
#' @param path File path location of CFinder.
#' 
#' @return  GeneModules = n x 3 dimensional data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
#' 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
findModules.CFinder.once <- function(adj, path, min.module.size, i){
  
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'undirected', weighted = T, diag = F)
  
  # Get modules using CFinder
  system(paste0('rm -rf ./tmp',i))
  system(paste0('mkdir ./tmp',i))
  
  ## Write network as an edgelist to a file
  elist = igraph::as_edgelist(g)
  elist = cbind(elist, 1)
  file = paste0('./tmp',i,'/input.txt')
  utils::write.table(elist, file = file, row.names = F, col.names = F, quote=F, sep = '\t')
  
  ## Run CFinder
  system(paste(paste0(path, 'CFinder_commandline64'),
               '-l',paste0(path, 'licence.txt'),
               '-i',file))
  
  ## Get output 
  d = list.dirs(paste0('./tmp',i,'/input.txt_files/'), full.names = T)
  
  mod = utils::read.table(paste0(d[2],'/communities'), skip = 7, sep = '\n') 
  mod$moduleNumber = 1:dim(mod)[1]
  geneModules = plyr::ddply(mod, .variables = "moduleNumber", .fun = function(x){
    mod = data.frame(Gene.ID = stringr::str_split(x$V1, ':')[[1]][2] %>%
                       stringr::str_split(' ') %>% 
                       unlist %>%
                       unique %>%
                       setdiff(c('')))
    mod$moduleNumber = x$moduleNumber
    mod$moduleSize = length(unique(mod$Gene.ID))
    return(mod)
  })
  system(paste0('rm -rf ./tmp',i))
  
  # Get individual cluster assignment for each gene from the community object
  geneModules = geneModules %>%
    dplyr::group_by(.data$Gene.ID) %>%
    dplyr::top_n(1, .data$moduleSize) %>%
    dplyr::top_n(1, .data$moduleNumber) %>%
    dplyr::select(-.data$moduleSize) %>%
    dplyr::mutate(moduleNumber = factor(.data$moduleNumber),
                  moduleNumber = as.numeric(.data$moduleNumber))
  
  # Add missing genes
  Gene.ID = setdiff(igraph::V(g)$name, geneModules$Gene.ID)
  geneModules = rbind(data.frame(geneModules), 
                      data.frame(Gene.ID = Gene.ID,
                                 moduleNumber = 0))
  
  # Rename modules with size less than min module size to 0
  filteredModules = geneModules %>% 
    dplyr::group_by(.data$moduleNumber) %>%
    dplyr::summarise(counts = length(unique(.data$Gene.ID))) %>%
    dplyr::filter(.data$counts >= min.module.size)
  geneModules$moduleNumber[!(geneModules$moduleNumber %in% filteredModules$moduleNumber)] = 0
  
  # Change cluster number to color labels
  geneModules$moduleNumber = as.numeric(factor(geneModules$moduleNumber))
  geneModules$moduleLabel = WGCNA::labels2colors(geneModules$moduleNumber)
  
  return(unique(geneModules))
}