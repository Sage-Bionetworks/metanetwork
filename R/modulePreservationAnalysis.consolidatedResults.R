#### Code to consolidate module conservation analysis ####
# reference: cell type markers
# test: rank consensus network across different sparsity and sparrow2Bonferroni networks from different methods

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

## Needs the dev branch
library(rGithubClient)

# Login to synapse
synapseLogin()

# Github links
thisRepo <- getRepo(repository = "th1vairam/metanetwork", ref="branch", refName="moduleAnal")
thisFile1 <- getPermlink(repository = thisRepo, repositoryPath = 'R/modulePreservationAnalysis.SGE.R')
thisFile2 <- getPermlink(repository = thisRepo, repositoryPath = 'R/modulePreservationAnalysis.cellMarkers.R')
thisFile3 <- getPermlink(repository = thisRepo, repositoryPath = 'R/modulePreservationAnalysis.consolidatedResults.R')

# Write reulst to synapse folder
parentId = "syn4974056"

# Synapse store metadata
activityName = "Module Preservation Analysis"
activityDescription = "Module Preservation Analysis against Cell type specific markers as sub-networks"

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
collectGarbage()

# Generate submission scripts for each comaprison
for (name in rownames(Data.Files)){
  # Get folder name to retrive results
  folderName = paste(Data.Files[name,c('file.method','file.sparsityMethod','file.disease')], collapse='_')
  folderName = paste(getwd(), 
                     paste('cellMarkers','as_ref',folderName,'as_test',sep='_'), 
                     sep='/')	

  # Read main results
  Main = fread(paste0(folderName, '/Main.tsv'), data.table=F, header = T)

  # Read random permutation results
  Rand = lapply(1:100, function(i, folderName){ try({results = fread(paste0(folderName, '/Rand.', i, '.tsv'),data.table=F,header=T)}, silent = T)}, folderName)

  # Remove failed runs
  Rand[which(sapply(Rand, length) == 1)] = NULL

  # Combine all permutation runs
  Rand = plyr::join_all(Rand, by = "moduleName")

  # Calculate z-scores
  propertyName = c("modTest", "cor.Adj", "meanAdj", "cor.PCor", "meanPCor", "changePCor", "meankIM", "changekIM2kALL")
  zstats = sapply(propertyName, function(prop, Rand, Main){
    z = (Main[, grep(prop, colnames(Main))] - rowMeans(Rand[, grep(prop, colnames(Rand))])) / apply(Rand[, grep(prop, colnames(Rand))], 1, sd)
  }, Rand, Main)
  colnames(zstats) = paste('Z', colnames(zstats), sep = '.') 
  zstats = zstats %>% data.frame 
  zstats$Z.min = apply(dplyr::select(zstats, Z.modTest, Z.cor.Adj, Z.meanAdj, Z.cor.PCor, Z.changekIM2kALL), 1, min)
  zstats$Z.median = apply(dplyr::select(zstats, Z.modTest, Z.cor.Adj, Z.meanAdj, Z.cor.PCor, Z.changekIM2kALL), 1, median) 
  zstats = cbind(Main[,1:2], zstats)

  # Synapse metadata
annotateList = list(dataType  = 'metaData',
                    sparsityMethod = Data.Files[name,'file.sparsityMethod'],
                    method = Data.Files[name,'file.method'],
                    disease = Data.Files[name,'file.disease'],
                    tissueType  = 'DLPFC',
                    normalization	= 'None',
                    testNet	= paste(Data.Files[name,c('file.method','file.sparsityMethod','file.disease')], collapse='_'), 
                    refNet	= 'CellMarkers',
                    organism	= 'HomoSapiens',
                    moduleMethod	= 'igraph:fast_greedy')
used = c(Data.Files[name,'file.id'], Module.Files[name, 'file.id'], 'syn4259377', 'syn4259379', 'syn4893059')

# Create results folder
foldObj = Folder(name = paste(Data.Files[name,c('file.method','file.sparsityMethod','file.disease')], collapse='_'), parentId = parentId)
annotations(foldObj) = annotateList
foldObj = synStore(foldObj)

# Create results metric file
modPresMetric = File(paste0(folderName, '/Main.tsv'), name = "Module Preservation Metrics", parentId = foldObj$properties$id)
annotations(modPresMetric) = annotateList
annotations(modPresMetric)$fileType = "tsv"
modPresMetric = synStore(modPresMetric, activityName = activityName, activityDescription = activityDescription, used = used, executed = list(thisFile1, thisFile2, thisFile3))
                     
# Create permuted results metric file
write.table(Rand, file = paste0(folderName, '/Rand.tsv'), sep = '\t', row.names=F, quote=F)
modPermPresMetric = File(paste0(folderName, '/Rand.tsv'), name = "Permuted Module Preservation Metrics", parentId = foldObj$properties$id)
annotations(modPermPresMetric) = annotateList
annotations(modPresMetric)$fileType = "tsv"
modPermPresMetric = synStore(modPermPresMetric, activityName = activityName, activityDescription = activityDescription, used = used, executed = list(thisFile1, thisFile2, thisFile3))

# Create zstats file
write.table(zstats, file = paste0(folderName, '/zstats.tsv'), sep = '\t', row.names=F, quote=F)
zstatObj = File(paste0(folderName, '/zstats.tsv'), name = "Module Preservation Zscores", parentId = foldObj$properties$id)
annotations(zstatObj) = annotateList
annotations(modPresMetric)$fileType = "tsv"
zstatObj = synStore(zstatObj, activityName = activityName, activityDescription = activityDescription, used = used, executed = list(thisFile1, thisFile2, thisFile3))

}
