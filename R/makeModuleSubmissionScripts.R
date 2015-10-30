#!usr/bin/env Rscript

# Submission Script in R
# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/home/ec2-user/Work/Github/metanetwork/R')

# Load libraries
library(plyr)
library(dplyr)
library(synapseClient)

# login to synapse
synapseLogin()

# Get all files and folder
Data.Files = synQuery('select * from file where 
                      projectId=="syn2397881" and 
                      fileType == "rda" and 
                      method == "rankconsensus" and 
                      sparsityMethod != "correlationBonferroni" and
                      sparsityMethod != "correlationFDR" and 
                      sparsityMethod != "wgcna"', blockSize = 100)
Data.Files = Data.Files$collectAll()
Data.Files = Data.Files %>%
  dplyr::mutate(uniqueName = paste(tools::file_path_sans_ext(file.name), file.disease, file.tissueType, sep = "."))

Module.Files = synQuery('select * from file where 
                        projectId=="syn2397881" and 
                        fileType == "tsv" and 
                        method == "rankconsensus" and
                        moduleMethod == "igraph:fast_greedy" and
                        sparsityMethod != "correlationBonferroni" and
                        sparsityMethod != "correlationFDR" and 
                        sparsityMethod != "wgcna"', blockSize = 100)
Module.Files = Module.Files$collectAll()
Module.Files = Module.Files %>%
  dplyr::filter(is.na(file.enrichmentMethod), file.name %in% grep('fast_greedy Modules', Module.Files$file.name, value=T)) %>%
  tidyr::separate(file.name, into = c("file.name","method","analysis"), sep = " ") %>%
  dplyr::mutate(uniqueName = paste(file.name, file.disease, file.tissueType, sep = "."))

All.Files = filter(Data.Files, !(uniqueName %in% Module.Files$uniqueName))

# Make directory and write shell scripts for running these files
system('mkdir sgeModuleSubmissions')
fp_all = file(paste('sgeModuleSubmissions/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)
for (id in All.Files$file.id){
  fp = file (paste('/home/ec2-user/Work/Github/metanetwork/R/sgeModuleSubmissions/SUB',id,sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /home/ec2-user/Work/Github/metanetwork/R/getModules.fastGreedy.R',id), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('sgeModuleSubmissions/allSubmissions.sh'),'a+')    
  cat(paste('qsub','-cwd','-V','-l h_vmem=7G', '-l mem_free=7G', paste('/home/ec2-user/Work/Github/metanetwork/R/sgeModuleSubmissions/SUB',id,sep='.'),
            '-o',paste('/home/ec2-user/Work/Github/metanetwork/R/sgeModuleSubmissions/SUB',id,'o',sep='.'),
            '-e',paste('/home/ec2-user/Work/Github/metanetwork/R/sgeModuleSubmissions/SUB',id,'e',sep='.')),
      file = fp_all,
      sep = '\n')
  close(fp_all)
}

