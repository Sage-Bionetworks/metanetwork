## Quick hack in sge to calculate module preservation stats with randomisation

# Load require libraries
library(synapseClient)
library(data.table)
library(igraph)
library(WGCNA)
library(plyr)
library(dplyr)

tmp = available.packages()
if (any(tmp[,1] %in% "glasso")){
	library('glasso')
} else {
	install.packages('glasso'); 
	library('glasso')
}

synapseLogin()

# Read command line arguments
args = commandArgs(TRUE)

# Load input test RData
load(args[1])

# Get folder names to save
folderName = args[2]

# Get filenames to save 
fileName = args[3]

# Load input reference RData from synapse
load(synGet('syn4973033')@filePath)

# Set working directory
setwd(folderName)

#############################################################################################
#### Function to calculate partial correlation coefficients ####
partialCorrelation <- function(C, adjMatrix){
  
#   browser()
  ind = arrayInd(which(adjMatrix == 0), dim(adjMatrix))
  gl = glasso(C, rho = 0, zero = ind)
  
  n = ncol(adjMatrix)
  pc = -gl$wi/sqrt(matrix(rep(diag(gl$wi), n),n,n)*t(matrix(rep(diag(gl$wi), n),n,n)))
  colnames(pc) = colnames(adjMatrix)
  rownames(pc) = rownames(adjMatrix)
  
#   print('Completed pc successfully')
  return(pc)
}

#### Function to calculate module preservation metrics ####
modulePreservationMetrics <- function(x, refModLabels, testModLabels, refNet, testNet, refExp, testExp){
  
  print(x)
  refModuleNames = names(refModLabels)[refModLabels == x]
  moduleSize = sum(refModLabels == x)
  
  ## Extract sub network modules 
  sgRef = igraph::subgraph(refNet, refModuleNames)
  sgTest = igraph::subgraph(testNet, refModuleNames)
  
  ## Extract expression matrix
  refExp = refExp[refExp$ensembl_gene_id %in% igraph::V(sgRef)$name,-(2)]
  rownames(refExp) = refExp$ensembl_gene_id
  refExp = refExp[igraph::V(sgRef)$name,-(1)]
  
  testExp = testExp[testExp$ensembl_gene_id %in% igraph::V(sgTest)$name,-(2)]
  rownames(testExp) = testExp$ensembl_gene_id
  testExp = testExp[igraph::V(sgTest)$name,-(1)]
  
  ## Calculate modularity of test network
  modTest = igraph::modularity(sgTest, unclass(factor(testModLabels[refModuleNames])))
  
  ### Density preservation stats
  # Correlation between adjacency matrices
  adjRef = igraph::get.adjacency(sgRef)
  adjTest = igraph::get.adjacency(sgTest)
  adjTest = adjTest[rownames(adjRef), colnames(adjRef)]
  cor.Adj = WGCNA::bicor(as.vector(adjRef), as.vector(adjTest),  use = 'pairwise.complete.obs')
  
  # Module density in test network
  meanAdj = igraph::graph.density(sgTest)
    
  ### Correlation preservation stats
  # Get correlation matrix
  refCor = WGCNA::bicor(t(refExp),  use = 'pairwise.complete.obs')
  testCor = WGCNA::bicor(t(testExp),  use = 'pairwise.complete.obs')
  
  try({
  # Get partial correlation matrix
  refCor = partialCorrelation(refCor, adjRef)
  testCor = partialCorrelation(testCor, adjTest)
  }, silent = T)
  
  # Correlation between partial correlation matrices
  cor.PCor = WGCNA::bicor(as.vector(refCor), as.vector(testCor), use = 'pairwise.complete.obs')
  
  # Mean partial correlation in test network
  meanPCor = mean(abs(testCor))
  
  # Fold change in partial correlation
  changePCor = meanPCor/mean(abs(refCor))
  
  ### Connectivity preservation stats
  # Degree of nodes (k.all)
  testDegreeAll = igraph::degree(testNet)
  
  # Degree of all nodes in sub-netowrk (k.in)
  testDegree = igraph::degree(sgTest)
  
  # Calculate mean intra modular connectivity
  meankIM = mean(testDegree)
  changekIM2kALL = mean(testDegree/testDegreeAll[names(testDegree)], na.rm=T)
  
  return(data.frame(
    moduleSize = moduleSize,
    modTest = modTest,
    cor.Adj = cor.Adj,
    meanAdj = meanAdj,
    cor.PCor = cor.PCor, 
    meanPCor = meanPCor,
    changePCor = changePCor, 
    meankIM = meankIM,
    changekIM2kALL = changekIM2kALL))
}
#############################################################################################

#############################################################################################
#### Data download and preparation ####
refNet = igraph::graph.adjacency(netData.ref$refNet, mode = 'undirected', weighted = NULL, diag = F)
testNet = igraph::graph.adjacency(netData.test$testNet, mode = 'undirected', weighted = NULL, diag = F)

refModLabels = netData.ref$refModLabels$modulelabels
names(refModLabels) = netData.ref$refModLabels$GeneIDs

# Assign modules wiht less than 20 genes to NoModule
testModLabels = netData.test$testModLabels
testModLabels = split(testModLabels, testModLabels$modulelabels)
tmp = sapply(testModLabels, dim)
testModLabels[which(tmp[1,]<20)] = lapply(testModLabels[which(tmp[1,]<20)], function(x){
  x$moduleNumber = 0
  x$modulelabels = "NoModule"
  return(x)
})
tmp = ldply(testModLabels)[,-(1)]
testModLabels = tmp$modulelabels
names(testModLabels) = tmp$GeneIDs

refExp = netData.ref$refExp
testExp = netData.test$testExp

# Perform randomisation (if needed)
if (fileName != 'Main'){
  print('Performing randomization...')
  ind = sample(vcount(refNet))
  adjRefNet = get.adjacency(refNet)
  refModLabels = refModLabels[rownames(adjRefNet)]
  rownames(adjRefNet)[ind] = rownames(adjRefNet)
  colnames(adjRefNet)[ind] = colnames(adjRefNet)
  permRefNet = igraph::graph.adjacency(adjRefNet, mode = 'undirected', weighted = NULL, diag = F)
  
  permModLabels = refModLabels
  names(permModLabels)[ind] = names(refModLabels)
  collectGarbage()
  
  ind = sample(vcount(testNet))
  adjTestNet = get.adjacency(testNet)
  testModLabels = testModLabels[rownames(adjTestNet)]
  rownames(adjTestNet)[ind] = rownames(adjTestNet)
  colnames(adjTestNet)[ind] = colnames(adjTestNet)
  permTestNet = igraph::graph.adjacency(adjTestNet, mode = 'undirected', weighted = NULL, diag = F)
  
  # Package actual data and submit them to sge
  refNet = permRefNet;
  testNet = permTestNet;
  refModLabels = permModLabels                                    
}

refModules = setdiff(unique(refModLabels), 'NoModule')
testModules = setdiff(unique(testModLabels), 'NoModule')

presMetrics = lapply(refModules, modulePreservationMetrics, refModLabels, testModLabels, refNet, testNet, refExp, testExp)

presMetrics = plyr::ldply(presMetrics)
presMetrics = cbind(data.frame(moduleName = refModules), presMetrics)
#############################################################################################

#############################################################################################
#### Write results to file ####
write.table(presMetrics, file = paste(folderName, paste(fileName,'tsv',sep='.'),sep='/'), sep='\t', row.names=F, quote=F)

# Write completed file list
write.table(paste(folderName, paste(fileName,'tsv',sep='.'),sep='/'), file = paste(folderName, 'Completed.txt', sep='/'), sep='\n', append=T, quote=F, col.names=F, row.names=F)
#############################################################################################
