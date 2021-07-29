
# Function to get modules from network adjacency matrix using MEGENA algorithm
findModules.megena <- function(adj, method="pearson", FDR = 0.05, module.pval = 0.05, hub.pval=0.05){
  
  library(MEGENA)
  method = method # method for correlation. either pearson or spearman. 
  FDR.cutoff = FDR # FDR threshold to define significant correlations upon shuffling samples. 
  module.pval = module.pval # module significance p-value. Recommended is 0.05. 
  hub.pval = hub.pval # connectivity significance p-value based random tetrahedral networks
  cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
  hub.perm = 100; # number of permutations for calculating connectivity significance p-value. 
  
  