modulePreservationMetrics <- function(x, refModLabels, refNet, testNet){
  
  moduleSize = sum(refModLabels == x)
  
  sgRef = induced_subgraph(refNet, names(refModLabels)[refModLabels == x])
  sgTest = induced_subgraph(testNet, names(refModLabels)[refModLabels == x])
  
  ## Module properties in test network
  meanAdj = graph.density(sgTest)
  meanClCoeff = transitivity(sgTest, type = "global")
  
  ## Connectivity preservation stats
  # Correlation between adjacency matrices
  adjRef = as_adjacency_matrix(sgRef)
  adjTest = as_adjacency_matrix(sgTest)
  adjTest = adjTest[rownames(adjRef), colnames(adjRef)]
  cor.Adj = cor(as.vector(adjRef), as.vector(adjTest))
  
  # Correlation between intramodular degree
  refDegree = degree(sgRef)
  testDegree = degree(sgTest)
  cor.kIM = bicor(refDegree, testDegree[names(refDegree)])
  
  # Correlation between clustering coefficients
  transRef = transitivity(sgRef, v=V(sgRef), type = 'barrat')
  names(transRef) = V(sgRef)$name
  transTest = transitivity(sgTest, v=V(sgTest), type = 'barrat')
  names(transTest) = V(sgTest)$name
  cor.Cl.Coef = cor(transRef, transTest[names(transRef)], use='p')
  
  ## Expression preservation stats
  
  
  return(data.frame(
    moduleSize = moduleSize,
    meanAdj = meanAdj,
    meanClCoeff = meanClCoeff,
    cor.Adj = cor.Adj,
    cor.kIM = cor.kIM,
    cor.Cl.Coef = cor.Cl.Coef))
}