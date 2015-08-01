library(synapseClient)
library(data.table)
library(rGithubClient)

synapseLogin()

# Filter gene sets for background
GL_OBJ = synGet('syn4868362')
load(GL_OBJ@filePath)

# Get background gene list
EXP_OBJ = synGet('syn4259377')
EXP = fread(EXP_OBJ@filePath, data.table=F)
backGroundGenes = unique(EXP$V2[-(1)])

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
                       x = x[len>minSize && len<maxSize]
                       return(x)
                     },
                     minSize,
                     maxSize)
  len = sapply(GeneLists, length)
  GeneLists = GeneLists[len != 0]
  
  return(GeneLists)
}

# Filter gene sets
GeneSets = filterGeneSets(GeneSets, backGroundGenes)

thisFileName <- 'filterGeneSets.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/metanetwork", 
                    ref="branch", 
                    refName='enrich')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('R/', thisFileName))

# Write to synapse
save(list = 'GeneSets', file = 'GeneLists.RData')
OBJ <- File('./GeneLists.RData',
            name = 'Enrichr Merged Gene Sets in RList Format (filtered)',
            parentId = 'syn4597301')
OBJ <- synStore(OBJ, 
                used = c('syn4868362', 'syn4259377'), 
                activityName = 'Gene Sets Curation',
                executed = thisFile)
