# Function to get modules from network adjacency matrix using GANXiS community detection algorithm v3.0.2
findModules.GANXiS <- function(adj, path, nperm = 10, min.module.size = 30){
  
  # Note: For this function to work get the source software from syn7806859, unzip and supply the path for GANXiSw.jar
  
  # Input
  #      adj = n x n upper triangular adjacency in the matrix class format
  #      nperm = number of permutation on the gene ordering 
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
  
  # Make adjacency matrix symmetric
  adj = adj + t(adj)
  adj[diag(adj)] = 0
  
  # Compute modules by permuting the labels nperm times
  all.modules = plyr::llply(1:nperm, .fun= function(i, adj, path, min.module.size){
    # Permute gene ordering
    ind = sample(1:dim(adj)[1], dim(adj)[1], replace = FALSE)
    adj1 = adj[ind,ind]
    
    mod = NA; Q = NA; Qds = NA;
    tryCatch({
      # Find modules 
      mod = findModules.GANXiS.once(adj1, path, min.module.size)
    
      # Compute local and global modularity
      adj1[lower.tri(adj1)] = 0
      Q = compute.Modularity(adj1, mod)
      Qds = compute.ModularityDensity(adj1, mod)
    }, error = function(e){
      mod = NA; Q = NA; Qds = NA;
    })
    
    return(list(mod = mod, Q = Q, Qds = Qds))
  }, adj, path, min.module.size)
  
  # Find the best module based on Q and Qds
  tmp = plyr::ldply(all.modules, function(x){
    data.frame(Q = x$Q, Qds = x$Qds)
  }) %>%
    na.omit %>%
    dplyr::mutate(r = base::rank(Q, na.last = FALSE)+base::rank(Qds, na.last = FALSE))
  ind = which.max(tmp$r)
  
  mod = all.modules[[ind]]$mod
  
  return(mod)
}

findModules.GANXiS.once <- function(adj, path, min.module.size){
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'undirected', weighted = T, diag = F)
  
  # Get modules using GANXiS
  system('rm -rf ./tmp')
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