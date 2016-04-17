#!/usr/bin/env Rscript

# Function to perform enrichment analysis od modules (from synapse as RData file)
# Get arguments from comman line
args = commandArgs(TRUE)

# Clear R console screen output
cat("\014")

# Set library and working directories
.libPaths('/shared/rlibs')
setwd('/shared/Github/metanetwork/R')
############################################################################################################

############################################################################################################
#### Libraries ####
library(synapseClient)
library(dplyr)
library(WGCNA)
library(tools)
library(stringr)
library(igraph)
library(data.table)
library(biomaRt)
library(CovariateAnalysis)

# Needs the dev branch
library(githubr)

# Login to synapse                   
key = read.table('/shared/synapseAPIToken')
synapseLogin(username = 'th_vairam', apiKey = key$V1)
############################################################################################################

############################################################################################################
#### Github commit ####
# Get github links for provenance
thisFileName1 = 'enrichModules.R'
thisFileName2 = 'makeEnrichSubmissionScripts.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/metanetwork", 
                    ref="branch", 
                    refName='enrich')

thisFile1 <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('R/', thisFileName1))
thisFile2 <- getPermlink(repository = thisRepo,
                         repositoryPath=paste0('R/', thisFileName2))

# Synapse specific parameters
activityName = 'Module enrichment'
activityDescription = 'Enrichment analysis of network modules using Fishers exact test'
############################################################################################################

############################################################################################################
#### Function definitions ####
# Function to filter Gene Sets
filterGeneSets <- function(GeneLists, # List of lists
                           genesInBackground, # background set of genes
                           minSize = 10,
                           maxSize = 1000){
  GeneLists = lapply(GeneLists, 
                     function(x, genesInBackground){
                       x = lapply(x, 
                                  function(x, genesInBackground){
                                    return(intersect(x, genesInBackground))
                                  },
                                  genesInBackground)
                       return(x)
                     }, 
                     genesInBackground)
  
  GeneLists = lapply(GeneLists, 
                     function(x, minSize, maxSize){
                       len = sapply(x, length)
                       x = x[len>minSize & len<maxSize]
                       return(x)
                     },
                     minSize,
                     maxSize)
  len = sapply(GeneLists, length)
  GeneLists = GeneLists[len != 0]
  
  return(GeneLists)
}

# Function to perform Fishers enrichment analysis
fisherEnrichment <- function(genesInSignificantSet, # A character vector of differentially expressed or some significant genes to test
                             genesInGeneSet, # A character vector of genes in gene set like GO annotations, pathways etc...
                             genesInBackground # Background genes that are 
){
  genesInSignificantSet = intersect(genesInSignificantSet, genesInBackground) # back ground filtering
  genesInNonSignificantSet = base::setdiff(genesInBackground, genesInSignificantSet)
  genesInGeneSet = intersect(genesInGeneSet, genesInBackground) # back ground filtering
  genesOutGeneSet = base::setdiff(genesInBackground,genesInGeneSet)
  
  pval = fisher.test(
    matrix(c(length(intersect(genesInGeneSet, genesInSignificantSet)),             
             length(intersect(genesInGeneSet, genesInNonSignificantSet)),
             length(intersect(genesOutGeneSet, genesInSignificantSet)),
             length(intersect(genesOutGeneSet, genesInNonSignificantSet))), 
           nrow=2, ncol=2),
    alternative="greater")
  OR = (length(intersect(genesInGeneSet, genesInSignificantSet)) * length(intersect(genesOutGeneSet, genesInNonSignificantSet))) / (length(intersect(genesInGeneSet, genesInNonSignificantSet)) * length(intersect(genesOutGeneSet, genesInSignificantSet)))
  return(data.frame(pval = pval$p.value,
                    ngenes = length(genesInGeneSet),
                    noverlap = length(intersect(genesInGeneSet, genesInSignificantSet)),
                    Odds.Ratio = OR,
                    Genes = paste(intersect(genesInGeneSet, genesInSignificantSet), collapse = '|')
  )
  )
}
############################################################################################################

############################################################################################################
#### Get gene sets ####
# Download enrichr gene sets from synapse
GL_OBJ = synGet('syn5923958')
ALL_USED_IDs = GL_OBJ$properties$id
load(GL_OBJ@filePath)
############################################################################################################

############################################################################################################
#### Get modules ####
# Download modules from synapse
MOD_OBJ = synapseClient::synGet(args[1])
ALL_USED_IDs = c(ALL_USED_IDs, MOD_OBJ$properties$id)
FNAME = args[2]
parentId = synGet(MOD_OBJ$properties$parentId)$properties$parentId

# Load modules 
MOD = data.table::fread(MOD_OBJ@filePath, data.table=F, header=T)
############################################################################################################

############################################################################################################
#### Background gene list ####
# Convert ensemble gene id's to hgnc symbols using biomart
# ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# ensg2hgnc = getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', values = MOD$EnsembleID, mart = ensembl)
ensg2hgnc = downloadFile('syn5923981')
ALL_USED_IDs = c(ALL_USED_IDs, 'syn5923981')
backGroundGenes = unique(ensg2hgnc$hgnc_symbol)

MOD = merge(MOD, ensg2hgnc, by.x = 'EnsembleID', by.y = 'ensembl_gene_id', all.x=T)
############################################################################################################

############################################################################################################
#### Filter gene list ####
GeneSets = filterGeneSets(GeneSets, backGroundGenes, minSize = 10, maxSize = 5000)
############################################################################################################

############################################################################################################
#### Perform enrichment analysis ####
# Perform enrichment analysis (for modules greater than 20 genes only)
# Perform enrichment analysis (for modules greater than 20 genes only)
enrichResults = list()
for (name in unique(MOD$moduleLabel)){
  genesInModule = unique(MOD$hgnc_symbol[MOD$moduleLabel == name])  
  if (length(genesInModule) > 20){
    enrichResults[[name]] = lapply(GeneSets,
                                   function(x, genesInModule, genesInBackground){
                                     tmp = as.data.frame(t(sapply(x, fisherEnrichment, genesInModule, genesInBackground)))
                                     tmp = rownameToFirstColumn(tmp,'GeneSetName')
                                     return(tmp)
                                   },
                                   unique(genesInModule), unique(backGroundGenes)) %>%
      rbindlist(use.names=TRUE, idcol = 'Category') %>%
      dplyr::mutate(fdr = p.adjust(pval, 'fdr'))
  } else {
    enrichResults[[name]] = data.frame(GeneSetName = NA, pval = NA, ngenes = NA, noverlap = NA, OR = NA, Category = NA, fdr = NA)
  }
  writeLines(paste0('Completed ',name))  
}
gc()

# Write results to file
tmp = rbindlist(enrichResults, use.names = TRUE, idcol = 'ComparisonName', fill = TRUE) %>%
  filter(!is.na(Category))

tmp$pval = unlist(tmp$pval)
tmp$ngenes = unlist(tmp$ngenes)
tmp$noverlap = unlist(tmp$noverlap)
tmp$Odds.Ratio = unlist(tmp$Odds.Ratio)
tmp$Genes = unlist(tmp$Genes)

write.table(tmp, file = paste(gsub(' ','_',FNAME),'enrichmentResults.tsv',sep='_'), sep='\t', row.names=F)
gc()
############################################################################################################

############################################################################################################
#### Write to synapse ####
# Create a folder for module results
obj = Folder(name = 'Module Enrichment', parentId = parentId)
obj = synStore(obj)

# Write results to synapse
algo = 'Fisher'
ENR_OBJ = File(paste(gsub(' ','_',FNAME),'enrichmentResults.tsv',sep='_'), name = 'BIC Rank Consensus', parentId = obj$properties$id)
annotations(ENR_OBJ) = annotations(MOD_OBJ)
ENR_OBJ@annotations$fileType = 'tsv'
ENR_OBJ@annotations$analysisType = 'moduleEnrichment'
ENR_OBJ@annotations$enrichmentMethod = 'Fisher'
ENR_OBJ@annotations$enrichmentGeneSet = 'Enrichr, AD, SCZ, cranio, genefamily, cellMarkers'
ENR_OBJ@annotations$moduleMethod = NULL
ENR_OBJ@annotations$modularity = NULL
ENR_OBJ = synStore(ENR_OBJ, 
                   executed = list(thisFile1, thisFile2),
                   used = ALL_USED_IDs,
                   activityName = activityName,
                   activityDescription = activityDescription)
############################################################################################################

############################################################################################################
# Write completed files to synapse
write.table(ENR_OBJ$properties$id, file = 'CompletedEnrichmentIDs.txt', sep='\n', append=T, quote=F, col.names=F, row.names=F)
writeLines(paste('Completed',FNAME,'and stored in',ENR_OBJ$properties$id))
############################################################################################################
