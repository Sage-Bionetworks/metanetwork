#' Compute Hypergeometric PValue
#' 
#' Computes a hypergeometric enrichment test of geneset `geneSetToTest` in a 
#' given module of genes `moduleGenes` given a total background set of genes 
#' `backgroundGenes`. Returns a P-Value of the enrichment test.
#' 
#' @param geneSetToTest Required. A character vector of gene IDs to test for an 
#' enrichment in `moduleGenes`.
#' @param moduleGenes Required. A character vector of gene IDs to test for being 
#' enriched with `geneSetToTest`.
#' @param backgroundGenes Required. A character vector of gene IDs to serves as 
#' the background gene set.
#' 
#' @export
#' @return numerical PValue
#' 
calcHyperPval <- function(geneSetToTest, moduleGenes, backgroundGenes){
  q = length(intersect(geneSetToTest,moduleGenes))
  m = length(geneSetToTest)
  n = length(backgroundGenes) - length(geneSetToTest)
  k = length(moduleGenes)  
  p.val = stats::phyper(q,m,n,k,lower.tail = T)
  
  return(p.val)
}