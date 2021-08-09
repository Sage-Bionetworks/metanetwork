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
  cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
  hub.perm = 100; # number of permutations for calculating connectivity significance p-value. 
  col_adj <- colnames(adj)
  new_col_adj <- c()
  for (i in col_adj){
    temp <- as.character(i)
    temp <- strsplit(temp,'[|]')[[1]][1]
    new_col_adj <- c(new_col_adj,temp)
  }
  colnames(adj) <- new_col_adj
  rownames(adj) <- new_col_adj
  g  <- graph.adjacency(adj,weighted=TRUE)
  g <- get.data.frame(g)
  colnames(g) <- c('row','col','weight')
  run.par = doPar & (getDoParWorkers() == 1) 
  if (run.par){
    cl <- parallel::makeCluster(n.cores)
    registerDoParallel(cl)
    # check how many workers are there
    cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
    }
  el <- calculate.PFN(g,doPar = TRUE,num.cores = n.cores,
                      keep.track = FALSE)

  MEGENA.output <- do.MEGENA(g,
                             mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = FALSE,
                             min.size = 10,max.size = vcount(g)/2,
                             doPar = TRUE,num.cores = nc,n.perm = hub.perm,
                             save.output = FALSE)
  summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                         mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                         min.size = 10,max.size = vcount(g)/2,
                                         annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                         output.sig = TRUE)
  
  