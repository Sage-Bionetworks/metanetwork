#!/usr/bin/env Rscript

# Function to perform enrichment analysis od modules (from synapse as RData file)
# Get arguments from comman line
args = commandArgs(TRUE)

# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/home/ec2-user/Work/Github/metanetwork/R')
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

# Needs the dev branch
library(rGithubClient)

# login to synapse
synapseLogin()
############################################################################################################

############################################################################################################
#### Github commit ####
# Get github links for provenance
thisFileName = 'enrichModules.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/metanetwork", 
                    ref="branch", 
                    refName='enrich')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('R/', thisFileName))

# Synapse specific parameters
activityName = 'Module enrichment'
activityDescription = 'Enrichment analysis of network modules using hypergeometric method'
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
  return(data.frame(pval = pval$p.value,
                    ngenes = length(genesInGeneSet),
                    noverlap = length(intersect(genesInGeneSet, genesInSignificantSet))))
}

# Function to convert rownames to first column of a df
rownameToFirstColumn <- function(DF,colname){
  DF <- as.data.frame(DF)
  DF[,colname] <- row.names(DF)
  DF <- DF[,c(dim(DF)[2],1:(dim(DF)[2]-1))]
  return(DF)
}
############################################################################################################

############################################################################################################
#### Get gene sets ####
# Download enrichr gene sets from synapse
GL_OBJ = synGet('syn4867851')
ALL_USED_IDs = GL_OBJ$properties$id
load(GL_OBJ@filePath)

gsets = c("Achilles_fitness_decrease", "Achilles_fitness_increase", "Allen_Brain_Atlas_down", "Allen_Brain_Atlas_up",
          "BioCarta", "CMAP_down", "CMAP_up", "Cancer_Cell_Line_Encyclopedia", "ChEA", "Cross_Species_Phenotype",
          "Disease_Signatures_from_GEO_down", "Disease_Signatures_from_GEO_up", "Drug_Perturbations_from_GEO",
          "ENCODE_Histone_Modifications_2013", "ESCAPE", "GO_Biological_Process", "GO_Cellular_Component", "GO_Molecular_Function",
          "GeneSigDB", "Genome_Browser_PWMs.1", "HMDB_Metabolites", "HomoloGene", "Human_Gene_Atlas", "KEGG_2015",
          "MGI_Mammalian_Phenotype", "MGI_Mammalian_Phenotype_Level_3", "MGI_Mammalian_Phenotype_Level_4", "MSigDB_Computational",
          "MSigDB_Oncogenic_Signatures", "Mouse_Gene_Atlas", "NCI-60_Cancer_Cell_Lines", "NCI-Nature", 
          "NURSA_Human_Endogenous_Complexome", "OMIM_Disease", "OMIM_Expanded", "PPI_Hub_Proteins", "Pfam_InterPro_Domains",
          "Phosphatase_Substrates_from_DEPOD", "Reactome", "SILAC_Phosphoproteomics", "TF-LOF_Expression_from_GEO", 
          "TargetScan_microRNA", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB",
          "Transcription_Factor_PPIs", "Virus_Perturbations_from_GEO_down", "Virus_Perturbations_from_GEO_up", "WikiPathways_2015",
          "GeneFamily","CellMarkers")
GeneSets = GeneSets[gsets]
############################################################################################################

############################################################################################################
#### Get modules ####
# Download modules from synapse
MOD_OBJ = synapseClient::synGet(args[1])
ALL_USED_IDs = c(ALL_USED_IDs, MOD_OBJ$properties$id)
FNAME = str_replace_all(tools::file_path_sans_ext(MOD_OBJ$properties$name), 'Modules', 'Enrichment')
parentId = MOD_OBJ$properties$parentId

# Load modules 
MOD = data.table::fread(MOD_OBJ@filePath, data.table=F, header=T)
############################################################################################################

############################################################################################################
#### Background gene list ####
# Convert ensemble gene id's to hgnc symbols using biomart
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensg2hgnc = getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', values = MOD$GeneIDs, mart = ensembl)
backGroundGenes = unique(ensg2hgnc$hgnc_symbol)

MOD = merge(MOD, ensg2hgnc, by.x = 'GeneIDs', by.y = 'ensembl_gene_id', all.x=T)
############################################################################################################

############################################################################################################
#### Filter gene list ####
GeneSets = filterGeneSets(GeneSets, backGroundGenes)
############################################################################################################

############################################################################################################
#### Perform enrichment analysis ####
# Perform enrichment analysis (for modules greater than 20 genes only)
enrichResults = list()
for (name in unique(MOD$modulelabels)){
  genesInModule = unique(MOD$hgnc_symbol[MOD$modulelabels == name])  
  if (length(genesInModule) > 20){
    tmp = lapply(GeneSets,
                 function(x, genesInModule, genesInBackground){
                   tmp = as.data.frame(t(sapply(x, fisherEnrichment, genesInModule, genesInBackground)))
                   tmp$pval = p.adjust(tmp$pval,'fdr')
                   tmp = rownameToFirstColumn(tmp,'GeneSetName')
                   return(tmp)
                 },
                 unique(genesInModule), unique(backGroundGenes))
    
    for (name1 in names(tmp))
      tmp[[name1]]$category = name1
    
    enrichResults[[name]] = as.data.frame(rbindlist(tmp))
  } else {
    enrichResults[[name]] = data.frame(GeneSetName =NA, pval =NA, ngenes =NA, noverlap = NA, category = NA)
  }
  writeLines(paste0('Completed ',name))  
}

# Write results to file
for(name in names(enrichResults))
  enrichResults[[name]]$ComparisonName = name
enrichmentResults = as.data.frame(rbindlist(enrichResults))
enrichmentResults$ngenes = unlist(enrichmentResults$ngenes)
enrichmentResults$noverlap = unlist(enrichmentResults$noverlap)
enrichmentResults = enrichmentResults[enrichmentResults$pval<=0.5,]
write.table(enrichmentResults, file = paste(gsub(' ','_',FNAME),'enrichmentResults.tsv',sep='_'), sep='\t', row.names=F, quote=F)
collectGarbage()
############################################################################################################

############################################################################################################
#### Write to synapse ####
# Write results to synapse
algo = 'Fisher'
ENR_OBJ = File(paste(gsub(' ','_',FNAME),'enrichmentResults.tsv',sep='_'), name = paste(FNAME,algo), parentId = parentId)
annotations(ENR_OBJ) = annotations(MOD_OBJ)
ENR_OBJ@annotations$fileType = 'tsv'
ENR_OBJ@annotations$enrichmentMethod = 'Fisher'
ENR_OBJ = synStore(ENR_OBJ, 
                   executed = thisFile,
                   used = ALL_USED_IDs,
                   activityName = activityName,
                   activityDescription = activityDescription)
############################################################################################################

############################################################################################################
# Write completed files to synapse
write.table(ENR_OBJ$properties$id, file = 'CompletedEnrichmentIDs.txt', sep='\n', append=T, quote=F, col.names=F, row.names=F)
writeLines(paste('Completed',FNAME,'and stored in',ENR_OBJ$properties$id))
############################################################################################################
