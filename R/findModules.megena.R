# Function to get modules from network adjacency matrix
findModules.megena <- function(data, method = "pearson", FDR.cutoff = 0.05, module.pval = 0.05, hub.pval = 0.05,doPar = TRUE){
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
n.cores <- detectCores() - 1; # number of cores/threads to call for PCP
    doPar <-TRUE; # do we want to parallelize?
    method = method# method for correlation. either pearson or spearman. 
    FDR.cutoff = FDR.cutoff# FDR threshold to define significant correlations upon shuffling samples. 
    module.pval = module.pval # module significance p-value. Recommended is 0.05. 
    hub.pval = hub.pval # connectivity significance p-value based random tetrahedral networks
    cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
    hub.perm = 100; # number of permutations for calculating connectivity significance p-value. 
                                    
    ijw <- calculate.correlation(data,doPerm = cor.perm,output.corTable = FALSE,output.permFDR = FALSE)

    cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = "")
    el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores,keep.track = FALSE)
    g <- graph.data.frame(el,directed = FALSE)

    MEGENA.output <- do.MEGENA(g,
    mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
    min.size = 10,max.size = vcount(g)/2,
    doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
    save.output = FALSE)

    summary.output <- MEGENA.ModuleSummary(MEGENA.output,
    mod.pvalue = module.pval,hub.pvalue = hub.pval,
    min.size = 10,max.size = vcount(g)/2,
    annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
    output.sig = TRUE)


    megena_Res = megena_Res %>% filter(module.size > min_module - 1)
    gene_modules = data.frame('Gene.ID'= character(),'moduleNumber'=numeric(),'moduleSize'=numeric())


    for (i in 1:nrow(megena_Res)){
    
    genes =  megena_Res$module.hub[i]
    genes = as.character(genes)
    genes = strsplit(genes,',')[[1]]
    for (k in 1:length(genes)){
        genes[k] = strsplit(genes[k],"[(]")[[1]][1]
    }
    gene_modules = rbind(gene_modules,data.frame('Gene.ID'= genes,'moduleNumber'=as.character(rep(megena_Res$module.id[i],length(genes))),
                                        'moduleSize'= megena_Res$module.size[i]))
    }


    gene_modules = gene_modules %>%
    group_by(Gene.ID) %>%
    dplyr::top_n(1, moduleSize) %>%
    dplyr::top_n(1, moduleNumber) %>%
    dplyr::select(-moduleSize) %>%
    dplyr::mutate(moduleNumber = factor(moduleNumber),
                    moduleNumber = as.numeric(moduleNumber))

    filteredModules = gene_modules %>% 
    dplyr::group_by(moduleNumber) %>%
    dplyr::summarise(counts = length(unique(Gene.ID))) %>%
    dplyr::filter(counts >= 30)
    gene_modules$moduleNumber[!(gene_modules$moduleNumber %in% filteredModules$moduleNumber)] = 0



    # Change cluster number to color labels
    gene_modules$moduleNumber = as.numeric(factor(gene_modules$moduleNumber))
    gene_modules$moduleLabel = WGCNA::labels2colors(gene_modules$moduleNumber)
    mod = gene_modules[c('Gene.ID','moduleNumber','moduleLabel')]
    return(mod)

}