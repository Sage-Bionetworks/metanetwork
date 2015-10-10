#### Code to delploy module preservation analysis ####
# reference: AD, MCI, NCI - Saprrow2Bonferroni
# test: AD, MCI, NCI - Saprrow2Bonferroni

## It is assumed your working directory is where this file is

setwd('/home/ec2-user/Work/Github/metanetwork/R')

# Clear R console screen output
cat("\014")  

# Load require libraries
library(synapseClient)
library(data.table)
library(igraph)
library(org.Hs.eg.db)
library(annotate)
library(tools)
library(biomaRt)
library(rGithubClient)
library(WGCNA)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(knitr)
library(stringr)

# Login to synapse
synapseLogin()

####################################################################################################################################
#### Function definitions ####
filterModuleLabels <-  function(refModLabels, n = 20){
  mod2remove = refModLabels %>% 
    group_by(modulelabels) %>% 
    summarise(count = n()) %>% 
    dplyr::filter(count < 20) %>% 
    dplyr::select(modulelabels) %>%
    unlist
  refModLabels$moduleNumber[refModLabels$modulelabels %in% mod2remove] = 0
  refModLabels$modulelabels[refModLabels$modulelabels %in% mod2remove] = 'NoModule'
  return(refModLabels)
}
####################################################################################################################################

#### Download expression data ####
EXP_ID = 'syn4259377'
EXP = fread(synGet(EXP_ID)@filePath, data.table = F, header=T)

#### Download expression data ####
COV_ID = 'syn4259379'
COVARIATES = fread(synGet(COV_ID)@filePath, data.table = F, header = T)
COVARIATES = split(COVARIATES, COVARIATES$cogdx)

#### Separate expression based on cogdx ####
exp = lapply(COVARIATES, function(x,EXP) { x = EXP[, colnames(EXP) %in% x$SampleID]; x = cbind(EXP[,'ensembl_gene_id', drop=F],x)}, EXP)
exp = exp[c(4,2,1)]
names(exp) = c("AD","MCI","NCI")
collectGarbage()

#### Get reference network: cell type markers ####
# Download AD related gene sets from synapse
GL_OBJ = synGet('syn4893059');
load(GL_OBJ@filePath) # this will load a list named GeneSets

# Get human related mapping
Hs = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org") # use this one when biomart.org is down
Hs = useDataset("hsapiens_gene_ensembl", Hs)
NGeneSets = lapply(GeneSets[grep('Zhang',names(GeneSets),value=T)], function(x, Hs){  
  human_ensg2symbol = getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                            filters = "hgnc_symbol",                         
                            values = x,
                            mart = Hs)
  return(human_ensg2symbol[,'ensembl_gene_id'])
}, Hs)

# Download adjacency matrices from synapse
Data.Files = synQuery('select * from file where projectId=="syn2397881" and fileType == "rda" and method == "rankconsensus" and sparsityMethod == "sparrow2Bonferroni"')
rownames(Data.Files) = Data.Files$file.disease

# Download modules from synapse
Module.Files = synQuery('select * from file where projectId=="syn2397881" and fileType == "tsv" and method == "rankconsensus" and moduleMethod == "igraph:fast_greedy" and sparsityMethod == "sparrow2Bonferroni"')
Module.Files = unique(Module.Files[grep('Modules',Module.Files$file.name),])
rownames(Module.Files) = Module.Files$file.disease
collectGarbage()

# Generate submission scripts for each comaprison
for (name.i in rownames(Data.Files)){
  for (name.j in rownames(Data.Files)){
    if (name.i != name.j){
      # Create folder to save files
      folderName = paste(getwd(), 
                         paste(name.i,'as_ref',name.j,'as_test',sep='_'), 
                         sep='/')
      system(paste('mkdir',folderName))
      
      # Download reference adjacency matrix and formulate an igraph object
      load(synGet(Data.Files[name.i,'file.id'])@filePath)
      
      # Download reference modules
      refModLabels = fread(synGet(Module.Files[name.i,'file.id'])@filePath, data.table=F, header=T)
      refModLabels = filterModuleLabels(refModLabels, n = 20)
      
      # Get reference expression data
      refExp = exp[[Data.Files[name.i,'file.disease']]]
      
      # Package actual data and submit them to sge
      netData = list(refNet = sparseNetwork, refModLabels = refModLabels, refExp = refExp)                                    
      collectGarbage()
      
      # Download test adjacency matrix and formulate an igraph object
      load(synGet(Data.Files[name.j,'file.id'])@filePath)
      
      # Download test modules
      testModLabels = fread(synGet(Module.Files[name.j,'file.id'])@filePath, data.table=F, header=T)
      testModLabels = filterModuleLabels(testModLabels, n = 20)
      
      # Get test expression data
      testExp = exp[[Data.Files[name.j,'file.disease']]]
      
      # Package actual data and submit them to sge
      netData = c(netData, list(testNet = sparseNetwork, testModLabels = testModLabels, testExp = testExp))                                    
      save(list = 'netData', file = paste(folderName, 'Input.RData',sep='/'))
      collectGarbage()
      
      # Track all subission scripts in one shell script
      fp_all = file(paste0(folderName, '/allSubmissions.sh'),'w')    
      cat('#!/bin/bash',file=fp_all,sep='\n')
      close(fp_all)
      
      # Create main submission script
      fp = file (paste0(folderName,'/Main.sh'), "w")
      cat('#!/bin/bash',
          'sleep 30',
          paste('Rscript','/home/ec2-user/Work/Github/metanetwork/R/modulePreservationAnalysis.SGE.R','Input.RData',folderName,'Main'),
          file = fp,
          sep = '\n')
      close(fp)
      
      # Add submission script to allSubmission list
      fp_all = file(paste0(folderName, '/allSubmissions.sh'),'a+')
      cat(paste('qsub','-cwd','-V',paste(folderName,'Main.sh',sep='/'),
                '-o',paste(folderName,'Main.o',sep='/'),
                '-e',paste(folderName,'Main.e',sep='/'),
                '-l mem=7GB'),
          file=fp_all,
          sep='\n')
      close(fp_all)
      
      # Create random networks for sge submission
      for (i in 1:150){
        # Create main submission script
        fp = file (paste(folderName, paste('Rand',i,'sh',sep='.'),sep='/'), "w")
        cat('#!/bin/bash',
            'sleep 30',
            paste('Rscript','/home/ec2-user/Work/Github/metanetwork/R/modulePreservationAnalysis.SGE.R',paste(folderName, 'Input.RData',sep='/'),folderName,paste('Rand',i,sep='.')),
            file = fp,
            sep = '\n')
        close(fp)
        
        # Add submission script to allSubmission list
        fp_all = file(paste0(folderName, '/allSubmissions.sh'),'a+')
        cat(paste('qsub','-cwd','-V',paste(folderName, paste('Rand',i,'sh',sep='.'),sep='/'),
                  '-o',paste(folderName, paste('Rand',i,'o',sep='.'),sep='/'),
                  '-e',paste(folderName, paste('Rand',i,'e',sep='.'),sep='/'),
                  '-l mem=7GB'),
            file=fp_all,
            sep='\n')
        close(fp_all)
      }
    }
  }
}
