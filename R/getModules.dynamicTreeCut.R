#!/usr/bin/env Rscript

# Function to get modules from network adjacency matrix (from synapse as rda file)
# Get arguments from comman line
args = commandArgs(TRUE)

# Function parameters to calculate modules
TOMType = 'unsigned' #other options are 'signed'
TOMDenom = 'min' #other options are 'mean',
linkageType = 'ward' #
distanceType = 'euclidean' #
minClusterSize = 30
deepSplit = F

# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/home/ec2-user/Work/Github/metanetwork/R')

# Load libraries
library(synapseClient)
library(dplyr)
library(WGCNA)
library(Rclusterpp)
library(tools)
library(stringr)

# Needs the dev branch
library(rGithubClient)

# Login to synapse
synapseLogin()

# Get github links for provenance
thisFileName = 'getModules.dynamicTreeCut.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/metanetwork", 
                    ref="branch", 
                    refName='modules')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('R/', thisFileName))

# Synapse specific parameters
activityName = 'Module Identification'
activityDescription = 'Clustering network genes in to modules using TOM and dynamic tree cut methodology'

# Get network from synapse (rda format)
NET_OBJ = synapseClient::synGet(args[1])
FNAME = tools::file_path_sans_ext(NET_OBJ$properties$name)
parentId = NET_OBJ$properties$parentId
  
# Load sparse network
load(NET_OBJ@filePath)
  
# Convert lsparseNetwork to adjacency
adjMat = as.matrix(sparseNetwork)*1

# Check adjMat is of class matrix
if (!is(adjMat, 'matrix'))
  stop('adjMat should be of class matrix')
  
# Get topological overlap matrix
TOM = WGCNA::TOMdist(adjMat, 
                     TOMType, 
                     TOMDenom, 
                     verbose = 2)
colnames(TOM) = colnames(adjMat)
rownames(TOM) = rownames(adjMat)
  
# Get hierarchical clusters from TOM
collectGarbage()
Rclusterpp::Rclusterpp.setThreads(threads=4)
clust = Rclusterpp::Rclusterpp.hclust(TOM,
                                      method = linkageType,
                                      distance = distanceType)
  
# Get individual clusters from the hierarchical tree
collectGarbage()
clust.numLabels = try(dynamicTreeCut::cutreeDynamic(clust,
                                                    minClusterSize = minClusterSize,
                                                    method = 'tree',
                                                    deepSplit = deepSplit))
  
if (is(clust.numLabels,'try-error'))
  clust.numLabels = rep(0, dim(adjMat)[1])
  
# Change cluster number to color labels
collectGarbage()
labels = WGCNA::labels2colors(clust.numLabels)
  
# Get results
geneModules = data.frame(GeneIDs = rownames(adjMat),
                         moduleNumber = clust.numLabels, 
                         modulelabels = labels,
                         row.names = rownames(adjMat))
  
MOD = list(geneModules = geneModules,
           TOM = TOM,
           clust = clust)

# Write results to synapse
save(list=c('MOD'), file = paste(FNAME,'RData',sep='.'))
MOD_OBJ = File(paste(FNAME,'RData',sep='.'), name = paste(FNAME,'TOM','hclust','Modules'), parentId = parentId)
annotations(MOD_OBJ) = annotations(NET_OBJ)
MOD_OBJ@annotations$moduleMethod = paste('WGCNA','dynamicTreeCut',sep=':')
MOD_OBJ@annotations$networkStorageType = 'full'
MOD_OBJ@annotations$fileType = 'RData'
MOD_OBJ@annotations$moduleParameters = 'linkage:ward; distance:eucledian; treecut:dynamictree; minsize:30; deepSplit:F'
MOD_OBJ1 = synStore(MOD_OBJ,
                   executed = thisFile, 
                   used = NET_OBJ,
                   activityName = activityName,
                   activityDescription = activityDescription)

write.table(MOD$geneModules, paste(FNAME,'tsv',sep='.'), sep='\t', row.names=F, quote=F)
MOD_OBJ = File(paste(FNAME,'tsv',sep='.'), name = paste(FNAME,'dynamicTreeCut','Modules'), parentId = parentId)
annotations(MOD_OBJ) = annotations(NET_OBJ)
MOD_OBJ@annotations$fileType = 'tsv'
MOD_OBJ@annotations$moduleParameters = 'linkage:ward; distance:eucledian; treecut:dynamictree; minsize:30; deepSplit:F'
MOD_OBJ@annotations$moduleMethod = paste('WGCNA','dynamicTreeCut',sep=':')
MOD_OBJ2 = synStore(MOD_OBJ, 
                   executed = thisFile,
                   used = NET_OBJ,
                   activityName = activityName,
                   activityDescription = activityDescription)

write.table(c(MOD_OBJ1$properties$id, MOD_OBJ2$properties$id), file = 'CompletedIDs.txt',append=T, quote=F, row.names=F, col.names=F)
writeLines(paste('Completed',FNAME,'and stored in',MOD_OBJ1$properties$id,MOD_OBJ2$properties$id))
