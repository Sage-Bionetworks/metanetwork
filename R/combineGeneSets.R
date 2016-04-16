#!/usr/bin/env Rscript

# Clear R console screen output
cat("\014")

# Set library and working directories
.libPaths('/shared/rlibs')
setwd('/shared/Github/metanetwork/R')

# Load libraries
library(synapseClient)
library(dplyr)
library(WGCNA)
library(tools)
library(stringr)
library(igraph)
library(data.table)
library(CovariateAnalysis)

# Needs the dev branch
library(githubr)

# Login to synapse                   
key = read.table('/shared/synapseAPIToken')
synapseLogin(username = 'th_vairam', apiKey = key$V1)

#### Github commit ####
# Get github links for provenance
thisFileName = 'combineGeneSets.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/metanetwork", 
                    ref="branch", 
                    refName='enrich')

thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('R/', thisFileName))

#### Get gene sets ####
# Download enrichr gene sets from synapse
GL_OBJ = synGet('syn4867851')
ALL_USED_IDs = GL_OBJ$properties$id
load(GL_OBJ@filePath) # This RData file will load a list named GeneSets

gsets = c("Allen_Brain_Atlas_down", "Allen_Brain_Atlas_up", "BioCarta_2015", "Cancer_Cell_Line_Encyclopedia", "ChEA", "Chromosome_Location", "Cross_Species_Phenotype",
          "Disease_Signatures_from_GEO_down", "Disease_Signatures_from_GEO_up", "Drug_Perturbations_from_GEO",
          "ENCODE_Histone_Modifications_2015", "ENCODE_TF_ChIP-seq_2015", "ESCAPE", "Epigenomics_Roadmap_HM_ChIP-seq",
          "GO_Biological_Process", "GO_Cellular_Component", "GO_Molecular_Function", "GeneSigDB", "Genome_Browser_PWMs.1",
          "HomoloGene", "HumanCyc", "Human_Gene_Atlas", "Human_Phenotype_Ontology", "KEA", "KEGG_2015", "MSigDB_Oncogenic_Signatures",
          "Mouse_Gene_Atlas", "OMIM_Disease","OMIM_Expanded", "PPI_Hub_Proteins", "Panther", "Reactome_2015", 
          "TF-LOF_Expression_from_GEO", "TRANSFAC_and_JASPAR_PWMs", "TargetScan_microRNA", "Transcription_Factor_PPIs","WikiPathways_2015")
GeneSets.enrichr = GeneSets[gsets]

# Download gene family genesets from synapse
allFiles <- synQuery('select * from file where parentId=="syn3240583"') %>%
  filter(file.name %in% c('GENEFAMILY'))
ALL_USED_IDs = c(ALL_USED_IDs, allFiles$file.id)
GeneSets.GeneFamily = downloadFile(allFiles$file.id) %>%
  dlply(.(Name), .fun = function(x){
    str_split(x$symBeforeOverlap, '\\|')[[1]]
  })

# Download AD related gene sets from synapse
GL_OBJ = synGet('syn4893059');
ALL_USED_IDs = c(ALL_USED_IDs, GL_OBJ$properties$id)
load(GL_OBJ@filePath)

GeneSets.Cell_Markers = GeneSets[c("Zhang:Astrocyte", "Zhang:Endothelial", "Zhang:Microglia", "Zhang:MyelinOligos",
                                   "Zhang:Neuron", "Zhang:NewOligos", "Zhang:OPC")]
GeneSets.AD = GeneSets[c("AD:GeneticLoci", "MouseMicroglia:2month_TG_vs_WT", "MouseMicroglia:4month_TG_vs_WT",
                         "MouseMicroglia:4month_vs_2month_TG-WT", "MouseMicroglia:6month_TG_vs_WT",
                         "MouseMicroglia:6month_vs_4month_TG-WT", "MouseMicroglia:8month_TG_vs_WT",
                         "MouseMicroglia:8month_vs_6month_TG-WT", "ROSMAP.NCIvsAD", "ROSMAP.MCIvsAD")]

# Download schizophrenia related gene sets from synapse
allFiles <- synQuery('select * from file where parentId=="syn3240583"') %>%
  filter(file.name %in% c('ASD.GENETICS', 'SCZ.GENETICS', 'GENEANNOT.ASD', 'GENEANNOT.SCZ'))
ALL_USED_IDs = c(ALL_USED_IDs, allFiles$file.id)
GeneSets.SCZ = sapply(allFiles$file.id, function(id){
  tmp = downloadFile(id) %>%
    dlply(.(Name), .fun = function(x){
      str_split(x$symBeforeOverlap, '\\|')[[1]]
    })
  tmp = tmp[setdiff(names(tmp), grep('Zhang', names(tmp), value = T))]
}) %>%
  do.call(c, .)

# Download cranio related gene sets from synapse
GL_OBJ = synGet('syn5752718')
ALL_USED_IDs = c(ALL_USED_IDs, GL_OBJ$properties$id)
load(GL_OBJ@filePath)

GeneSets.Cranio = GeneSets

# Combine all gene sets
GeneSets = c(GeneSets.enrichr, 
             list(Alzheimers = GeneSets.AD, 
                  Schizophrenia = GeneSets.SCZ, 
                  Cranio = GeneSets.Cranio, 
                  Genefamily = GeneSets.GeneFamily, 
                  Cell_Markers = GeneSets.Cell_Markers))

### Store results in synapse ###
save(list = 'GeneSets', file = 'metaNetGeneSets.RData')
obj = File('metaNetGeneSets.RData', name = 'Genesets for metanetworks', parentId = 'syn5923956')
obj = synStore(obj, used = ALL_USED_IDs, executed = thisFile, activity = 'Combine all related genesets')
