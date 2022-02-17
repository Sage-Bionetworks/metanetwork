#' Fishers Enrichment Analysis
#'
#'  This function to perform Fishers enrichment analysis. Tests the enrichment of
#'  `genesInGeneSet` in `genesInModule` considering a bacckground gene set of 
#'  `genesInBackground`.
#'
#' @param genesInModule Required. A character vector of differentially expressed genes or genes in a given module to test
#' @param genesInGeneSet Required. A character vector of genes in a gene set like GO annotations, pathways etc.
#' @param genesInBackground Required. A character vector of background genes
#' 
#' @return A dataframe consisting of column value of fisher test P-Value. number 
#' of genes in the test set, enrichment set, the overlap, odds ratio and the genes
#' from the test set which are present in the enrichment set.
#'
#' @export
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
  
  pval = stats::fisher.test(confusion.matrix, alternative="greater")
  OR = (confusion.matrix[1,1] * confusion.matrix[2,2])/
    (confusion.matrix[1,2] * confusion.matrix[2,1])
  
  return(data.frame(pval = pval$p.value,
                    nTested = length(genesInModule),
                    nInSet = length(genesInGeneSet),
                    noverlap = confusion.matrix[1,1],
                    OR = OR,
                    Genes = paste(intersect(genesInGeneSet, genesInSignificantSet), collapse = '|')) )
}