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
Module.Files = synQuery('select * from file where projectId=="syn2397881" and fileType == "tsv" and moduleMethod == "igraph:fast_greedy"', blockSize = 100)
Module.Files = Module.Files$collectAll()
Module.Files = Module.Files[is.na(Module.Files$file.enrichmentMethod),]
#Enrich.Files = synQuery('select * from file where projectId=="syn2397881" and fileType == "tsv" and enrichmentMethod == "Fisher" and enrichmentGeneSet == "AD"', blocksize = 100)
#Enrich.Files = Enrich.Files$collectAll()

# Make directory and write shell scripts for running these files
system('mkdir sgeEnrichSub')
fp_all = file(paste('./sgeEnrichSub/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)
for (id in Module.Files$file.id){
  fp = file (paste('/home/ec2-user/Work/Github/metanetwork/R/sgeEnrichSub/SUB',id,sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /home/ec2-user/Work/Github/metanetwork/R/enrichModules.R',id), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('./sgeEnrichSub/allSubmissions.sh'),'a+')    
  cat(paste('qsub','-cwd','-V',paste('/home/ec2-user/Work/Github/metanetwork/R/sgeEnrichSub/SUB',id,sep='.'),
            '-o',paste('/home/ec2-user/Work/Github/metanetwork/R/sgeEnrichSub/SUB',id,'o',sep='.'),
            '-e',paste('/home/ec2-user/Work/Github/metanetwork/R/sgeEnrichSub/SUB',id,'e',sep='.')),
      file=fp_all,
      sep='\n')
  close(fp_all)
}
