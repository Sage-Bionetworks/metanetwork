## Quick hack in sge to calculate module preservation stats with randomisation

# Load require libraries
library(synapseClient)
library(data.table)
library(igraph)
library(WGCNA)
library(dplyr)

# Read command line arguments
args = commandArgs(TRUE)
load(args[1])
folderName = args[2]
fileName = args[3]

# Set working directory
setwd(folderName)

#############################################################################################
#### Function to calculate module preservation metrics ####
modulePreservationMetrics <- function(x, refModLabels, refNet, testNet, refExp, testExp){
  
  refModuleNames = names(refModLabels)[refModLabels == x]
  moduleSize = sum(refModLabels == x)
  
  ## Extract sub network modules 
  sgRef = igraph::subgraph(refNet, refModuleNames)
  sgTest = igraph::subgraph(testNet, refModuleNames)
  
  ## Extract expression matrix
  refExp = refExp[refExp$ensembl_gene_id %in% igraph::V(sgRef)$name,-(1)]
  testExp = testExp[testExp$ensembl_gene_id %in% igraph::V(sgTest)$name,-(1)]
  
  refCor = WGCNA::bicor(t(refExp))
  testCor = WGCNA::bicor(t(testExp))
  
  cor.Cor = WGCNA::bicor(as.vector(refCor), as.vector(testCor))
  
  ## Connectivity preservation stats
  # Correlation between adjacency matrices
  adjRef = igraph::as_adjacency_matrix(sgRef)
  adjTest = igraph::as_adjacency_matrix(sgTest)
  adjTest = adjTest[rownames(adjRef), colnames(adjRef)]
  cor.Adj = WGCNA::bicor(as.vector(adjRef), as.vector(adjTest))
  
  # Correlation between intramodular degree
  refDegree = igraph::degree(sgRef)
  testDegree = igraph::degree(sgTest)
  cor.kIM = WGCNA::bicor(refDegree, testDegree[names(refDegree)])
   
  ## Module properties in test network
  meanAdj = igraph::graph.density(sgTest)
  meankIM = mean(testDegree)
  
  ## Fold change in module properties
  changeAdj = meanAdj/igraph::graph.density(sgRef)
  
  return(data.frame(
    moduleSize = moduleSize,
    meanAdj = meanAdj,
    meankIM = meankIM,    
    changeAdj = changeAdj,  
    cor.Cor = cor.Cor,
    cor.Adj = cor.Adj,
    cor.kIM = cor.kIM))
}
#############################################################################################

#############################################################################################
#### Data download and preparation ####
refNet = netData$refNet
testNet = netData$testNet
refModLabels = netData$refModLabels
testModLabels = netData$testModLabels
refExp = netData$refExp
testExp = netData$testExp
  
refModules = setdiff(unique(refModLabels), 'NoModule')
testModules = setdiff(unique(testModLabels), 'NoModule')
  
presMetrics = lapply(refModules, modulePreservationMetrics, refModLabels, refNet, testNet, refExp, testExp)
  
presMetrics = plyr::ldply(presMetrics)
presMetrics = cbind(data.frame(moduleName = refModules), presMetrics)
#############################################################################################

#############################################################################################
#### Write results to file ####
write.table(presMetrics, file = paste(folderName, paste(fileName,'tsv',sep='.'),sep='/'), sep='\t', row.names=F, quote=F)

# Write completed file list
write.table(paste(folderName, paste(fileName,'tsv',sep='.'),sep='/'), file = paste(folderName, 'Completed.txt', sep='/'), sep='\n', append=T, quote=F, col.names=F, row.names=F)
#############################################################################################

