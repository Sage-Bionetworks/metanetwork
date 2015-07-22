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
library(tools)
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
All.Files = synQuery('select name,id,disease from file where projectId=="syn2397881" and fileType == "rda"')
Finished.Files = synQuery('select name,id,disease from file where projectId=="syn2397881" and fileType == "RData"')

All.Files = All.Files[!(paste(tools::file_path_sans_ext(All.Files$file.name),All.Files$file.disease) %in%
                          paste(sapply(Finished.Files$file.name, function(x){strsplit(x," ")[[1]][1]}), Finished.Files$file.disease)),]

# Synapse specific parameters
activityName = 'Module Identification'
activityDescription = 'Clustering network genes in to modules using TOM and dynamic tree cut methodology'
  
for(i in 1:length(All.Files$file.id)){
  NET_OBJ = synGet(All.Files$file.id[i])
  FNAME = tools::file_path_sans_ext(NET_OBJ$properties$name)
  parentId = NET_OBJ$properties$parentId
  
  # Load sparse network
  load(NET_OBJ@filePath)
  
  # Convert lsparseNetwork to adjacency
  NET = as.matrix(sparseNetwork)*1
  
  # Get modules
  MOD = getModules(NET)
  
  # Write results to synapse
  save(list=c('MOD'), file = paste(FNAME,'RData',sep='.'))
  MOD_OBJ = File(paste(FNAME,'RData',sep='.'), name = paste(FNAME,'TOM','hclust','Modules'), parentId = parentId)
  annotations(MOD_OBJ) = annotations(NET_OBJ)
  MOD_OBJ@annotations$networkStorageType = 'full'
  MOD_OBJ@annotations$fileType = 'RData'
  MOD_OBJ@annotations$moduleParameters = 'linkage:ward_distance:eucledian_treecut:dynamictree_minsize:30_deepSplit:F'
  MOD_OBJ = tryCatch(synStore(MOD_OBJ,
                              executed = list(thisFile,moduleFun), 
                              used = NET_OBJ,
                              activityName = activityName,
                              activityDescription = activityDescription),
                     error = function(MOD_OBJ){
                       MOD_OBJ = synStore(MOD_OBJ,
                                          executed = list(thisFile,moduleFun), 
                                          used = NET_OBJ,
                                          activityName = activityName,
                                          activityDescription = activityDescription)
                       return(MOD_OBJ)
                       }
                     )
  
  write.table(MOD$geneModules, paste(FNAME,'tsv',sep='.'), sep='\t', row.names=F, quote=F)
  MOD_OBJ = File(paste(FNAME,'tsv',sep='.'), name = paste(FNAME,'Modules'), parentId = parentId)
  annotations(MOD_OBJ) = annotations(NET_OBJ)
  MOD_OBJ@annotations$fileType = 'tsv'
  MOD_OBJ@annotations$moduleParameters = 'linkage:ward_distance:eucledian_treecut:dynamictree_minsize:30_deepSplit:F'
  MOD_OBJ = tryCatch(synStore(MOD_OBJ,
                              executed = list(thisFile,moduleFun), 
                              used = NET_OBJ,
                              activityName = activityName,
                              activityDescription = activityDescription),
                     error = function(MOD_OBJ){
                       MOD_OBJ = synStore(MOD_OBJ,
                                          executed = list(thisFile,moduleFun), 
                                          used = NET_OBJ,
                                          activityName = activityName,
                                          activityDescription = activityDescription)
                       return(MOD_OBJ)
                     }
  )                     
  writeLines(paste('Completed',FNAME,'and stored in',MOD_OBJ$properties$id))
}