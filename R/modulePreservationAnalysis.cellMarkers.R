#### Code to delploy module preservation analysis ####
# reference: cell type markers
# test: rank consensus network across different sparsity and sparrow2Bonferroni networks from different methods

## It is assumed your working directory is where this file is

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
names(exp) = c("AD","MCI","NC\I")

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
Data.Files1 = synQuery('select * from file where projectId=="syn2397881" and fileType == "rda" and method == "rankconsensus" and sparsityMethod != "correlationBonferroni" and sparsityMethod != "correlationFDR" and sparsityMethod != "wgcna"')
Data.Files2 = synQuery('select * from file where projectId=="syn2397881" and fileType == "rda" and sparsityMethod == "sparrow2Bonferroni" and method != "correlationBonferroni" and method != "correlationFDR" and sparsityMethod != "wgcna"')
Data.Files = unique(rbind(Data.Files1, Data.Files2))
rownames(Data.Files) = paste(Data.Files$file.method, Data.Files$file.sparsityMethod, Data.Files$file.disease, sep='.')

# Download modules from synapse
Module.Files1 = synQuery('select * from file where projectId=="syn2397881" and fileType == "tsv" and method == "rankconsensus" and moduleMethod == "igraph:fast_greedy" and sparsityMethod != "correlationBonferroni" and sparsityMethod != "correlationFDR"')
Module.Files2 = synQuery('select * from file where projectId=="syn2397881" and fileType == "tsv" and sparsityMethod == "sparrow2Bonferroni" and moduleMethod == "igraph:fast_greedy" and method != "correlationBonferroni" and method != "correlationFDR"')
Module.Files = rbind(Module.Files1, Module.Files2)
Module.Files = unique(Module.Files[grep('Modules',Module.Files$file.name),])
rownames(Module.Files) = paste(Module.Files$file.method, Module.Files$file.sparsityMethod, Module.Files$file.disease, sep='.')

# Generate submission scripts for each comaprison
for (name in rownames(Data.Files)){
  # Download test adjacency matrix and formulate an igraph object
  load(synGet(Data.Files[name,'file.id'])@filePath)
  testNet = igraph::graph.adjacency(sparseNetwork, mode = 'undirected', weighted = NULL, diag = F)
  
  # Download test modules
  testModLabels = fread(synGet(Module.Files[name,'file.id'])@filePath, data.table=F, header=T)
    
  # Get test expression data
  testExp = exp[[Data.Files[name,'file.disease']]]
  
  # Make reference adjacency matrix
  refNetwork = matrix(0, dim(sparseNetwork)[1], dim(sparseNetwork)[2])
  colnames(refNetwork) = colnames(sparseNetwork)
  rownames(refNetwork) = rownames(sparseNetwork)
  
  # Make reference modules
  refModLabels = testModLabels
  refModLabels$moduleNumber = 0
  refModLabels$modulelabels = "NoModule"
    
  # Assign cell specific network and modules
  for (i in names(NGeneSets)){
    refNetwork[rownames(refNetwork) %in% NGeneSets[[i]], colnames(refNetwork) %in% NGeneSets[[i]]] = 1
    refModLabels$moduleNumber[refModLabels$GeneIDs %in% NGeneSets[i]] = 0
    refModLabels$modulelabels[refModLabels$GeneIDs %in% NGeneSets[i]] = str_replace_all(i,':','.')
  }
  diag(refNetwork) = 0
  refNet = igraph::graph.adjacency(refNetwork, mode = 'undirected', weighted = NULL, diag = F)
  
  # Get test expression data
  refExp = exp[[Data.Files[name,'file.disease']]]
  
  # Create folder to save files
  folderName = paste(Data.Files[name,c('file.method','file.sparsityMethod','file.disease')], collapse='_')
  folderName = paste(getwd(), 
                     paste('cellMarkers','as_ref',folderName,'as_test',sep='_'), 
                     sep='/')
  system(paste('mkdir',folderName))
    
  # Package actual data and submit them to sge
  netData = list(refNet = refNet, testNet = testNet, 
                 refModLabels = refModLabels, testModLabels = testModLabels, 
                 refExp = refExp, testExp = testExp)                                    
  save(list = 'netData', file = paste(folderName, 'Input.RData',sep='/'))
    
  # Track all subission scripts in one shell script
  fp_all = file(paste0(folderName, '/allSubmissions.sh'),'w+')    
  cat('#!/bin/bash',file=fp_all,sep='\n')
  close(fp_all)
  
  # Create main submission script
  fp = file (paste0(folderName,'/Main.sh'), "w+")
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
  for (i in 1:2){
    # Create main submission script
    fp = file (paste(folderName, paste('Rand',i,'sh',sep='.'),sep='/'), "w+")
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