#!usr/bin/env Rscript

# Submission Script in R
# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/home/ec2-user/Work/Github/metanetwork/R')

# Load libraries
library(synapseClient)

# login to synapse
synapseLogin()

# Get all files and folder
Data.Files = synQuery('select * from file where projectId=="syn2397881" and fileType == "rda"', blockSize = 100)
Data.Files = Data.Files$collectAll()

Module.Files = synQuery('select * from file where projectId=="syn2397881" and fileType == "tsv" and moduleMethod == "igraph:fast_greedy"', blockSize = 100)
Module.Files = Module.Files$collectAll()
Module.Files = Module.Files[is.na(Module.Files$file.enrichmentMethod),]

All.Files = Data.Files[!(paste(tools::file_path_sans_ext(Data.Files$file.name),Data.Files$file.disease) %in%
                         paste(sapply(Module.Files$file.name, function(x){strsplit(x," ")[[1]][1]}), Module.Files$file.disease)),]

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
