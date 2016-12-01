# Function to get modules from network adjacency matrix using CFinder algorithm
findModules.CFinder <- function(adj, path, min.module.size = 20){
  
  # Note: For this function to work get the source software from syn7806853,
  # unzip and supply the path for CFinder executable
  
  # Input
  #      adj = n x n adjacency matrix in the ltCMatrix format
  #      path = location of CFinder
  
  # Output (list of following elements)
  #      modules = data frame of dimension n x 3 with columns named GeneIDs, moduleNumber and moduleLabel
  
  # Error functions
  if(class(adj) != "ltCMatrix")
    stop('Adjacency matrix should be of class ltCMatrix')
  
  if(dim(adj)[1] != dim(adj)[2])
    stop('Adjacency matrix should be symmetric')
  
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'upper', weighted = NULL, diag = F)
  
  # Get modules using GANXiS
  system('mkdir ./tmp')
  
  ## Write network as an edgelist to a file
  elist = igraph::as_edgelist(g)
  elist = cbind(elist, 1)
  write.table(elist, file = './tmp/input.txt', row.names = F, col.names = F, quote=F, sep = '\t')
  
  ## Run CFinder
  system(paste(paste0(path, 'CFinder_commandline64'),
               '-l',paste0(path, 'licence.txt'),
               '-i','./tmp/input.txt'))
  
  ## Get output 
  d = list.dirs('./tmp/input.txt_files/', full.names = T)
  
  mod = read.table(paste0(d[2],'/communities'), skip = 7, sep = '\n') 
  mod$moduleNumber = 1:dim(mod)[1]
  geneModules = plyr::ddply(mod, .(moduleNumber), .fun = function(x){
    mod = data.frame(Gene.ID = str_split(x$V1, ':')[[1]][2] %>%
                       str_split(' ') %>% 
                       unlist %>%
                       unique %>%
                       setdiff(c('')))
    mod$moduleNumber = x$moduleNumber
    return(mod)
  })
  system('rm -rf ./tmp')
  
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
  
  return(unique(geneModules))
}