#!/usr/bin/env Rscript

# Function to perform enrichment analysis od modules (from synapse as RData file)
# Get arguments from comman line
args = commandArgs(TRUE)

# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/home/ec2-user/Work/Github/metanetwork/R')

# Load libraries
library(synapseClient)
library(dplyr)
library(WGCNA)
library(tools)
library(stringr)
library(igraph)

# Needs the dev branch
library(rGithubClient)

# login to synapse
synapseLogin()

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

# Get background gene list
EXP_OBJ = synGet('syn4259377')
ALL_USED_IDs = EXP_OBJ$properties$id
EXP = data.table::fread(EXP_OBJ@filePath, data.table=F, header=T)
backGroundGenes = unique(EXP$V2[-(1)])

# Get filtered gene list from synapse
GL_OBJ = synGet('syn4871865')
ALL_USED_IDs = c(ALL_USED_IDs, GL_OBJ@properties$id)
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

# Get modules network from synapse (RData format)
MOD_OBJ = synapseClient::synGet(args[1])
FNAME = str_replace_all(tools::file_path_sans_ext(MOD_OBJ$properties$name), 'Modules', 'Enrichment')
parentId = MOD_OBJ$properties$parentId

# Load sparse network
MOD = data.table::fread(MOD_OBJ@filePath, data.table=F)

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

rownameToFirstColumn <- function(DF,colname){
  DF <- as.data.frame(DF)
  DF[,colname] <- row.names(DF)
  DF <- DF[,c(dim(DF)[2],1:(dim(DF)[2]-1))]
  return(DF)
}
enrichResults = list()
for (name in unique(MOD$modulelabels)){
  genesInModule = unique(EXP$hgnc_symbol[EXP$ensembl_gene_id %in% MOD$GeneIDs[MOD$modulelabels == name]])  
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
  
  enrichResults[[name]] = do.call(rbind,tmp)
  writeLines(paste0('Completed ',name))  
}

# Write results to file
for(name in names(enrichResults))
  enrichResults[[name]]$ComparisonName = name
enrichmentResults = do.call(rbind,enrichResults)
enrichmentResults$ngenes = unlist(enrichmentResults$ngenes)
enrichmentResults$noverlap = unlist(enrichmentResults$noverlap)
write.table(enrichmentResults, file = 'enrichmentResults.tsv', sep='\t', row.names=F, quote=F)
collectGarbage()

# Write results to synapse
algo = 'Fisher'
ENR_OBJ = File('enrichmentResults.tsv', name = paste(FNAME,algo), parentId = parentId)
annotations(ENR_OBJ) = annotations(MOD_OBJ)
ENR_OBJ@annotations$fileType = 'tsv'
ENR_OBJ@annotations$enrichmentMethod = 'Fisher'
ENR_OBJ = synStore(ENR_OBJ, 
                   executed = thisFile,
                   used = ALL_USED_IDs,
                   activityName = activityName,
                   activityDescription = activityDescription)

# Write completed files to synapse
write.table(ENR_OBJ$properties$id, file = 'CompletedModuleIDs.txt', sep='\n', append=T, quote=F, col.names=F, row.names=F)
writeLines(paste('Completed',FNAME,'and stored in',ENR_OBJ$properties$id))