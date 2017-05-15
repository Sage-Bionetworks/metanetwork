# Function to get modules from network adjacency matrix using CFinder algorithm
findModules.CFinder <- function(adj, path, nperm = 10, min.module.size = 30){
  
  # Note: For this function to work get the source software from syn7806853,
  # unzip and supply the path for CFinder executable
  
  # Input
  #      adj = n x n upper triangular adjacency in the matrix class format
  #      path = location of CFinder
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
      mod = findModules.CFinder.once(adj1, path, min.module.size)
    
      # Compute local and global modularity
      adj1[lower.tri(adj1)] = 0
      Q = compute.Modularity(adj1, mod)
      Qds = compute.ModularityDensity(adj1, mod)
    }, error = function(e){
      mod = NA;Q=NA;Qds=NA;
    })
    
      return(list(mod = mod, Q = Q, Qds = Qds))
  }, adj, path, min.module.size)
  
  # Find the best module based on Q and Qds
  tmp = plyr::ldply(all.modules, function(x){
    data.frame(Q = x$Q, Qds = x$Qds)
  }) %>%
    na.omit() %>%
    dplyr::mutate(r = base::rank(Q, na.last = F)+base::rank(Qds, na.last = F))
  ind = which.max(tmp$r)
  
  mod = all.modules[[ind]]$mod
  
  return(mod)
}

findModules.CFinder.once <- function(adj, path, min.module.size){
  
  # Convert lsparseNetwork to igraph graph object
  g = igraph::graph.adjacency(adj, mode = 'undirected', weighted = T, diag = F)
  
  # Get modules using CFinder
  system('rm -rf ./tmp')
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
    mod = data.frame(Gene.ID = stringr::str_split(x$V1, ':')[[1]][2] %>%
                       stringr::str_split(' ') %>% 
                       unlist %>%
                       unique %>%
                       setdiff(c('')))
    mod$moduleNumber = x$moduleNumber
    mod$moduleSize = length(unique(mod$Gene.ID))
    return(mod)
  })
  system('rm -rf ./tmp')
  
  # Get individual cluster assignment for each gene from the community object
  geneModules = geneModules %>%
    group_by(Gene.ID) %>%
    dplyr::top_n(1, moduleSize) %>%
    dplyr::top_n(1, moduleNumber) %>%
    dplyr::select(-moduleSize) %>%
    dplyr::mutate(moduleNumber = factor(moduleNumber),
                  moduleNumber = as.numeric(moduleNumber))
  
  # Add missing genes
  Gene.ID = setdiff(igraph::V(g)$name, geneModules$Gene.ID)
  geneModules = rbind(data.frame(geneModules), 
                      data.frame(Gene.ID = Gene.ID,
                                 moduleNumber = 0))
  
  # Rename modules with size less than min module size to 0
  filteredModules = geneModules %>% 
    dplyr::group_by(moduleNumber) %>%
    dplyr::summarise(counts = length(unique(Gene.ID))) %>%
    dplyr::filter(counts >= min.module.size)
  geneModules$moduleNumber[!(geneModules$moduleNumber %in% filteredModules$moduleNumber)] = 0
  
  # Change cluster number to color labels
  geneModules$moduleNumber = as.numeric(factor(geneModules$moduleNumber))
  geneModules$moduleLabel = WGCNA::labels2colors(geneModules$moduleNumber)
  
  return(unique(geneModules))
}