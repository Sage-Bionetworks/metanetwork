# Function to get modules from network adjacency matrix using GANXiS community detection algorithm v3.0.2
findModules.GANXiS <- function(adj, path, min.module.size = 20){
  
  # Note: For this function to work get the source software from syn7806859, unzip and supply the path for GANXiSw.jar
  
  # Input
  #      adj = n x n adjacency matrix in the ltCMatrix format
  #      path = location of GANXiS.jar
  
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
  write.table(elist, file = './tmp/input.txt', row.names = F, col.names = F, quote=F)
  
  ## Run GANXiS
  system(paste('java','-jar',paste0(path, 'GANXiSw.jar'),
               '-i','./tmp/input.txt',
               '-d', './tmp',
               '-Sym','1',
               '-ov','0',
               '-Onc','1',
               '-seed', '123456789'))
  
  ## Get output 
  mod = read.table('./tmp/SLPAw_input_run1_r0.5_v3_T100.icpm.node-com.txt')
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