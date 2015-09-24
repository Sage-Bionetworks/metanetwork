modulePreservationMetrics <- function(x, refModLabels, refNet, testNet, refExp, testExp){
  
  refModuleNames = names(refModLabels)[refModLabels == x]
  moduleSize = sum(refModLabels == x)
  
  ## Extract sub network modules 
  sgRef = induced_subgraph(refNet, refModuleNames)
  sgTest = induced_subgraph(testNet, refModuleNames)
  
  ## Extract expression data
  testExp = filter(testExp, ensembl_gene_id %in% refModuleNames)
  testExp = testExp[,-(1)]
  
  refExp = filter(refExp, ensembl_gene_id %in% refModuleNames)
  refExp = refExp[,-(1)]
  
  ## Connectivity preservation stats
  # Correlation between adjacency matrices
  adjRef = as_adjacency_matrix(sgRef)
  adjTest = as_adjacency_matrix(sgTest)
  adjTest = adjTest[rownames(adjRef), colnames(adjRef)]
  cor.Adj = bicor(as.vector(adjRef), as.vector(adjTest))
  
  # Correlation between intramodular degree
  refDegree = degree(sgRef)
  testDegree = degree(sgTest)
  cor.kIM = bicor(refDegree, testDegree[names(refDegree)])
  
  # Correlation between clustering coefficients
  transRef = transitivity(sgRef, v=V(sgRef), type = 'barrat')
  names(transRef) = V(sgRef)$name
  transTest = transitivity(sgTest, v=V(sgTest), type = 'barrat')
  names(transTest) = V(sgTest)$name
  cor.Cl.Coef = bicor(transRef, transTest[names(transRef)], use='p')
  
  ## Expression preservation stats  
  
  ## Module properties in test network
  meanAdj = graph.density(sgTest)
  meanClCoeff = transitivity(sgTest, type = "global")
  meankIM = mean(testDegree)
  
  return(data.frame(
    moduleSize = moduleSize,
    meanAdj = meanAdj,
    meanClCoeff = meanClCoeff,
    meankIM = meankIM,
    cor.Adj = cor.Adj,
    cor.kIM = cor.kIM,
    cor.Cl.Coef = cor.Cl.Coef))
}