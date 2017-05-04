# Function to perform Fishers enrichment analysis
fisherEnrichment <- function(genesInModule, # A character vector of differentially expressed genes or genes in a given module to test
                             genesInGeneSet, # A character vector of genes in a gene set like GO annotations, pathways etc...
                             genesInBackground # Background genes
){
  # Background filtering
  genesInSignificantSet = base::intersect(genesInModule, genesInBackground) # back ground filtering
  genesInGeneSet = base::intersect(genesInGeneSet, genesInBackground) # back ground filtering
  
  # Get non significant set
  genesInNonSignificantSet = base::setdiff(genesInBackground, genesInSignificantSet)
  genesOutGeneSet = base::setdiff(genesInBackground,genesInGeneSet)
  
  confusion.matrix = matrix(c(length(intersect(genesInGeneSet, genesInSignificantSet)),             
                              length(intersect(genesInGeneSet, genesInNonSignificantSet)),
                              length(intersect(genesOutGeneSet, genesInSignificantSet)),
                              length(intersect(genesOutGeneSet, genesInNonSignificantSet))),
                            nrow=2, ncol=2)
  
  pval = fisher.test(confusion.matrix, alternative="greater")
  OR = (confusion.matrix[1,1] * confusion.matrix[2,2])/
    (confusion.matrix[1,2] * confusion.matrix[2,1])
  
  return(data.frame(pval = pval$p.value,
                    nTested = length(genesInModule),
                    nInSet = length(genesInGeneSet),
                    noverlap = confusion.matrix[1,1],
                    OR = OR,
                    Genes = paste(intersect(genesInGeneSet, genesInSignificantSet), collapse = '|')) )
  }