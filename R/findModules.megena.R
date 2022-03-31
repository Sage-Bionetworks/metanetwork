#' Find Modules with Megena Clustering
#' 
#' This function finds modules from network adjacency
#' matrix using MEGENA. 
#'
#' @param data Required. An n x n upper triangular adjacency in the matrix class 
#' format.
#' @param method Optional. Method for correlation. either pearson or spearman. 
#' (Default = "pearson")
#' @param FDR.cutoff Optional. FDR threshold to define significant correlations 
#' upon shuffling samples. (Default = 0.05)
#' @param module.pval Optional. Module significance p-value. Recommended is 0.05. 
#' (Default = 0.05)
#' @param hub.pval Optional. Connectivity significance p-value based random 
#' tetrahedral networks. (Default = 0.05)
#' @param doPar Optional. If parallelization of clusters is allowed (Default =TRUE)
#' @param n.cores Optional. The number of cores/threads to call for PCP. If NULL, 
#' n.cores = detectCores() - 1. (Default = NULL)
#' @param cor.perm Optional. Number of permutations for calculating FDRs for all 
#' correlation pairs. (Default = 10)
#' @param hub.perm Optional. number of permutations for calculating connectivity 
#' significance p-value. (Default = 100)
#' @param min_module Optional. minimum number of nodes/genes allowed for filtering
#'  (Default = 30)
#'
#' @return  GeneModules = n x 3 dimensional data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
#' 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
findModules.megena <- function(data, method = "pearson", FDR.cutoff = 0.05, 
                               module.pval = 0.05, hub.pval = 0.05, doPar = TRUE,
                               n.cores = NULL, cor.perm = 10, hub.perm = 100, min_module = 30 ){
  if(is.null(n.cores)){
    n.cores <- parallel::detectCores() - 1
  }
  #library(MEGENA)
  #library(igraph)
                              
  ijw <- MEGENA::calculate.correlation(data,
                                       doPerm = cor.perm, 
                                       output.corTable = FALSE,
                                       output.permFDR = FALSE)

  cat(paste( "number of cores to use:", foreach::getDoParWorkers(), "\n", sep = ""))
  el = MEGENA::calculate.PFN(
    ijw[,1:3],
    doPar = doPar,
    num.cores = n.cores,
    keep.track = FALSE)

  g = igraph::graph.data.frame(el,directed = FALSE)

  MEGENA.output <- MEGENA::do.MEGENA(g,
    mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
    min.size = 10,max.size = igraph::vcount(g)/2,
    doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
    save.output = FALSE)

  summary.output <- MEGENA::MEGENA.ModuleSummary(MEGENA.output,
  mod.pvalue = module.pval,hub.pvalue = hub.pval,
  min.size = 10,max.size = igraph::vcount(g)/2,
  output.sig = TRUE)
  
  megena_Res <- summary.output$module.table
  megena_Res = megena_Res %>% dplyr::filter(.data$module.size > min_module - 1)
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
    dplyr::group_by(.data$Gene.ID) %>%
    dplyr::top_n(1, .data$moduleSize) %>%
    dplyr::top_n(1, .data$moduleNumber) %>%
    dplyr::select(-.data$moduleSize) %>%
    dplyr::mutate(moduleNumber = factor(.data$moduleNumber),
                    moduleNumber = as.numeric(.data$moduleNumber))

  filteredModules = gene_modules %>%  
    dplyr::group_by(.data$moduleNumber) %>% 
    dplyr::summarise(counts = length(unique(.data$Gene.ID))) %>% 
    dplyr::filter(.data$counts >= 30)
  
  gene_modules$moduleNumber[!(gene_modules$moduleNumber %in% filteredModules$moduleNumber)] = 0

  # Change cluster number to color labels
  gene_modules$moduleNumber = as.numeric(factor(gene_modules$moduleNumber))
  gene_modules$moduleLabel = WGCNA::labels2colors(gene_modules$moduleNumber)
  mod = gene_modules[c('Gene.ID','moduleNumber','moduleLabel')]
  return(mod)

}
