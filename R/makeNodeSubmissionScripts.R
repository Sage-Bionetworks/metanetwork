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
Data.Files = synQuery('select * from file where projectId=="syn2397881" and fileType == "rda"')
Node.Files = synQuery('select id, name, disease, density from file where projectId=="syn2397881" and fileType == "tsv"')
Node.Files = Node.Files[!is.na(Node.Files$file.density),]

All.Files = Data.Files[!(paste(tools::file_path_sans_ext(Data.Files$file.name),Data.Files$file.disease) %in%
                         paste(sapply(Node.Files$file.name, function(x){strsplit(x," ")[[1]][1]}), Node.Files$file.disease)),]

# Make directory and write shell scripts for running these files
system('mkdir sgeNodeSubmissions')
fp_all = file(paste('sgeNodeSubmissions/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)
for (id in All.Files$file.id){
  fp = file (paste('/home/ec2-user/Work/Github/metanetwork/R/sgeNodeSubmissions/SUB',id,sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /home/ec2-user/Work/Github/metanetwork/R/nodeNetworkProperties.R',id), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('sgeNodeSubmissions/allSubmissions.sh'),'a+')    
  cat(paste('qsub','-cwd','-V',paste('/home/ec2-user/Work/Github/metanetwork/R/sgeNodeSubmissions/SUB',id,sep='.'),
            '-o',paste('/home/ec2-user/Work/Github/metanetwork/R/sgeNodeSubmissions/SUB',id,'o',sep='.'),
            '-e',paste('/home/ec2-user/Work/Github/metanetwork/R/sgeNodeSubmissions/SUB',id,'e',sep='.')),
      file=fp_all,
      sep='\n')
  close(fp_all)
}
