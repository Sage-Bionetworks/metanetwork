# Find hypergeometric pvalue
calcHyperPval <- function(geneSetToTest,
                          moduleGenes,                             
                          backgroundGenes){
  q = length(intersect(geneSetToTest,moduleGenes))
  m = length(geneSetToTest)
  n = length(backgroundGeneSet) - length(geneSetToTest)
  k = length(moduleGenes)  
  p.val = phyper(q,m,n,k,lower.tail = T)
  
  return(p.val)
}