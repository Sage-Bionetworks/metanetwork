library(synapseClient)
library(data.table)
library(org.Hs.eg.db)
library(annotate)
library(tools)
library(biomaRt)

synapseLogin()

# Download data from synapse 
MG_Files = synQuery('select * from file where parentId=="syn4883033"')
ALL_USED_IDs = MG_Files$file.id

MG = lapply(MG_Files$file.id, function(id){fread(synGet(id)@filePath, data.table=F, header=T)})
names(MG) = file_path_sans_ext(MG_Files$file.name)

# Get mouse related mapping
Mm = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org") # use this one when biomart.org is down
Mm = useDataset("mmusculus_gene_ensembl", Mm)
human_homologues = getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"),
                         filters = "with_entrezgene",                         
                         values=T,
                         mart = Mm)

mouse_entrz2ensg = getBM(attributes = c("entrezgene","ensembl_gene_id"),
                         filters = "with_entrezgene",                         
                         values=T,
                         mart = Mm)

mouse_mapping = merge(mouse_entrz2ensg, human_homologues, by = 'ensembl_gene_id', all=T)

# Get human related mapping
Hs = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org") # use this one when biomart.org is down
Hs = useDataset("hsapiens_gene_ensembl", Hs)
human_ensg2symbol = getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                         filters = "ensembl_gene_id",                         
                         values = mouse_mapping$hsapiens_homolog_ensembl_gene,
                         mart = Hs)
setnames(human_ensg2symbol, 'ensembl_gene_id', 'hsapiens_homolog_ensembl_gene')
mouse_human_mapping = merge(mouse_mapping, human_ensg2symbol, by = 'hsapiens_homolog_ensembl_gene', all=T)

GeneSets_mouse = lapply(MG, function(x, mouse_human_mapping){  
  entrz_id = x$EntrezGene[(x$FDR <= 0.05 & abs(x[,colnames(x) %in% c('foldchange','FDR')]) >= 1.5)]  
  gs = unique(mouse_human_mapping$hgnc_symbol[mouse_human_mapping$entrezgene %in% entrz_id])
  gs = gs[!is.na(gs)]
  return(gs)
},mouse_human_mapping)
names(GeneSets_mouse) = gsub('mouse_microglia_','MouseMicroglia:',names(GeneSets_mouse))

tmp = fread(synGet('syn4597305')@filePath, data.table=F)
ALL_USED_IDs = c(ALL_USED_IDs, 'syn4597305')
GeneSets_CM = lapply(tmp$Name, function(x,tmp){
  strsplit(tmp[tmp$Name %in% x,'symBeforeOverlap'], '\\|')[[1]]
}, tmp)
names(GeneSets_CM) = tmp$Name

tmp = fread(synGet('syn4891675')@filePath, data.table=F, header=F)
ALL_USED_IDs = c(ALL_USED_IDs, 'syn4891675')
GeneSets_GeneticLoci = list("AD:GeneticLoci" = tmp$V1)

GeneSets = c(GeneSets_CM, GeneSets_GeneticLoci, GeneSets_mouse)

thisFileName <- './curateADGeneSets.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/metanetwork", 
                    ref="branch", 
                    refName='enrich')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('code/Rmd/', thisFileName))

# Push list to synapse (in RData format)
save(list = 'GeneSets', file = 'MergedGeneSets.RData')
GS <- File('./MergedGeneSets.RData', name = 'Merged Gene Sets TG.Mouse Cellmarkers GeneticLoci (in RData format)', parentId = "syn4597301")
GS <- synStore(GS, used = ALL_USED_IDs, activityName = 'Merging Gene Sets', executed = thisFile)