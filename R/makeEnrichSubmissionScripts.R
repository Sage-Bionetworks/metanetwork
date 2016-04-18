#!usr/bin/env Rscript

# Submission Script in R
# Clear R console screen output
cat("\014")

# Set library and working directories
.libPaths('/shared/rlibs')
setwd('/shared/Github/metanetwork/R')

# Load libraries
library(synapseClient)
library(dplyr)
library(tidyr)

# Login to synapse                   
key = read.table('/shared/synapseAPIToken')
synapseLogin(username = 'th_vairam', apiKey = key$V1)

# Get module files
# Get finished modules
module.files = synQuery('select * from file where projectId == "syn5584871" and analysisType == "statisticalNetworkReconstruction" and method == "bic" and organism == "HomoSapiens" and moduleMethod == "igraph:fast greedy"') %>%
  dplyr::mutate(uniqueName = paste(file.disease, file.tissueTypeAbrv, file.study, file.cogdx, sep = '.'))

# Get finished enrichment files
enrich.files = synQuery('select * from file where projectId == "syn5584871" and analysisType == "moduleEnrichment" and method == "bic" and organism == "HomoSapiens"') %>%
  dplyr::mutate(uniqueName = paste(file.disease, file.tissueTypeAbrv, file.study, file.cogdx, sep = '.'))

# module.files = module.files %>%
#   filter(!(uniqueName %in% enrich.files$uniqueName)) 

# Make directory and write shell scripts for running these files
system('mkdir sgeEnrichSub')
fp_all = file(paste('./sgeEnrichSub/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)

for (i in 1:dim(module.files)[1]){
  id = module.files$file.id[i]
  FNAME = paste(module.files$uniqueName[i],'enrichment', sep = '.')
  
  fp = file (paste('/shared/Github/metanetwork/R/sgeEnrichSub/SUB',id,sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /shared/Github/metanetwork/R/enrichModules.R',id,FNAME), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('./sgeEnrichSub/allSubmissions.sh'),'a+')    
  cat(paste('qsub','-cwd','-V',paste('/shared/Github/metanetwork/R/sgeEnrichSub/SUB',id,sep='.'),
            '-o',paste('/shared/Github/metanetwork/R/sgeEnrichSub/SUB',id,'o',sep='.'),
            '-e',paste('/shared/Github/metanetwork/R/sgeEnrichSub/SUB',id,'e',sep='.')),
      file=fp_all,
      sep='\n')
  close(fp_all)
}
