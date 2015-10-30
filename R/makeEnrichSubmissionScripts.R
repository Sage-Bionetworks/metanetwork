#!usr/bin/env Rscript

# Submission Script in R
# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/home/ec2-user/Work/Github/metanetwork/R')

# Load libraries
library(synapseClient)
library(dplyr)
library(tidyr)

# login to synapse
synapseLogin()

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

Enrich.Files = synQuery('select * from file where 
                        projectId=="syn2397881" and 
                        fileType == "tsv" and 
                        method == "rankconsensus" and
                        moduleMethod == "igraph:fast_greedy" and
                        sparsityMethod != "correlationBonferroni" and
                        sparsityMethod != "correlationFDR" and 
                        sparsityMethod != "wgcna"', blockSize = 100)
Enrich.Files = Enrich.Files$collectAll()
Enrich.Files = Enrich.Files %>%
  dplyr::filter(!is.na(file.enrichmentMethod)) %>%
  tidyr::separate(file.name, into = c("file.name","moduleMethod","analysis","enrichMethod"), sep = " ") %>%
  dplyr::mutate(uniqueName = paste(file.name, file.disease, file.tissueType, sep = "."))

Module.Files = filter(Module.Files, !(uniqueName %in% Enrich.Files$uniqueName))

# Make directory and write shell scripts for running these files
system('mkdir sgeEnrichSub')
fp_all = file(paste('./sgeEnrichSub/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)
for (id in UEnrich.Files$file.id){
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