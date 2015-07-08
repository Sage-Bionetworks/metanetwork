# Entry point for module detection
## It is assumed your working directory is where this file is

# Clear R console screen output
cat("\014")

# Clear R workspace
rm(list=ls())

# Load libraries
library(synapseClient)
library(dplyr)
library(WGCNA)
library(Rclusterpp)
library(tool)
library(stringr)

## Needs the dev branch
library(rGithubClient)

synapseLogin()

# Source module identification script to global env
source('./getModules.R')

# Get github links for provenance
thisFileName = 'detectModules_EntryPoint.R'
moduleFunName = 'getModules.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/metanetwork", 
                    ref="branch", 
                    refName='master')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('R/', thisFileName))

moduleFun <- getPermlink(repository = thisRepo,
                         repositoryPath=paste0('R/', moduleFunName))

# Get all files and folder
tmp = synQuery('select * from entity where parentId=="syn4593591"')
All.Files = tmp[tmp$entity.nodeType == 'file',]
All.Folder = tmp[tmp$entity.nodeType == 'folder',]

activityName = 'Module Identification'
  
activityDescription = 'Clustering network genes in to modules using TOM and dynamic tree cut methodology'
  
for(i in 1:length(All.Files$entity.id)){
  NET_OBJ = synGet(All.Files$entity.id[i])
  FNAME = tools::file_path_sans_ext(NET_OBJ$properties$name)
  parentId = NET_OBJ$properties$parentId
  
  # Load sparse network
  load(NET_OBJ@filePath)
  
  # Convert lsparseNetwork to adjacency
  NET = as.matrix(sparseNetwork)*1
  
  # Get modules
  MOD = getModules(NET[1:3000,1:3000])
  
  # Write results to synapse
  save(list=c('MOD'), file = paste(FNAME,'RData',sep='.'))
  MOD_OBJ = File(paste(FNAME,'RData',sep='.'), name = paste(FNAME,'RData',sep=' '), parentId = parentId)
  annotations(MOD_OBJ) = annotations(NET_OBJ)
  MOD_OBJ@annotations$networkStorageType = 'full'
  MOD_OBJ@annotations$fileType = 'RData'
  MOD_OBJ@annotations$moduleParameters = 'linkage:ward_distance:eucledian_treecut:dynamictree_minsize:30_deepSplit:F'
  MOD_OBJ = synStore(MOD_OBJ,
                     executed = c(thisFile,moduleFun), 
                     used = NET_OBJ,
                     activityName = activityName,
                     activityDescription = activityDescription)
  
  write.table(MOD$geneModules, paste(FNAME,'tsv',sep='.'), sep='\t', row.names=F, quote=F)
  MOD_OBJ = File(paste(FNAME,'tsv',sep='.'), name = paste(FNAME,'tsv',sep=' '), parentId = parentId)
  annotations(MOD_OBJ) = annotations(NET_OBJ)
  MOD_OBJ@annotations$fileType = 'tsv'
  MOD_OBJ@annotations$moduleParameters = 'linkage:ward_distance:eucledian_treecut:dynamictree_minsize:30_deepSplit:F'
  MOD_OBJ = synStore(MOD_OBJ,
                     executed = c(thisFile,moduleFun), 
                     used = NET_OBJ,
                     activityName = activityName,
                     activityDescription = activityDescription)
}
