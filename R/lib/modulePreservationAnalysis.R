modulePreservationAnalysis <- function(refNet, testNet, refModLabels, testModLabels, refExp, testExp, n= 100){
  
  refModules = unique(refModLabels); #setdiff(unique(refModLabels), 'NoModule')
  testModules = unique(testModLabels); #setdiff(unique(testModLabels), 'NoModule')
  
  actualPresMetrics = lapply(refModules, modulePreservationMetrics, refModLabels, refNet, testNet, refExp, testExp)
  
  actualPresMetrics = plyr::ldply(actualPresMetrics)
  actualPresMetrics = cbind(data.frame(moduleName = refModules), actualPresMetrics)
  
  tmp = do.call(cbind, lapply(actualPresMetrics[,-(1:2)], rank))
  colnames(tmp) = paste(colnames(tmp), 'rank',sep='.')
  actualPresMetrics.rank = cbind(actualPresMetrics[,1:2,drop=F], 
                                 as.data.frame(tmp), 
                                 data.frame(median.rank = rank(apply(tmp,1,median))))
  
  randomPresMetrics = lapply(1:n, function(i, refModLabels, testModLabels, refNet, testNet, refExp, testExp){  
    browser()
    ind = sample(vcount(refNet))
    adjRefNet = as_adj(refNet)
    refModLabels = refModLabels[rownames(adjRefNet)]
    rownames(adjRefNet)[ind] = rownames(adjRefNet)
    colnames(adjRefNet)[ind] = colnames(adjRefNet)
    permRefNet = igraph::graph.adjacency(adjRefNet, mode = 'undirected', weighted = NULL, diag = F)
    
    permModLabels = refModLabels
    names(permModLabels)[ind] = names(refModLabels)
    
    ind = sample(vcount(testNet))
    adjTestNet = as_adj(testNet)
    testModLabels = testModLabels[rownames(adjTestNet)]
    rownames(adjTestNet)[ind] = rownames(adjTestNet)
    colnames(adjTestNet)[ind] = colnames(adjTestNet)
    permTestNet = igraph::graph.adjacency(adjTestNet, mode = 'undirected', weighted = NULL, diag = F)
    
    refModules = unique(permModLabels)
    
    permPresMetrics = lapply(refModules, modulePreservationMetrics, permModLabels, permRefNet, permTestNet, refExp, testExp)
    permPresMetrics = plyr::ldply(permPresMetrics)
    permPresMetrics = cbind(data.frame(moduleName = refModules), permPresMetrics)  
  }, refModLabels, testModLabels[2:4], refNet, testNet, refExp, testExp)

  randomPresMetrics = plyr::join_all(randomPresMetrics, by = 'moduleName')  
  randomPresMetrics[is.na(randomPresMetrics)] = 0
  
  Zstatistics = lapply(c("meanAdj", "meanClCoeff", "cor.Adj", "cor.kIM", "cor.Cl.Coef"), function(x, actualPresMetrics, randomPresMetrics){
    ind = grep(x,colnames(randomPresMetrics))
    return(z = (actualPresMetrics[,x] - rowMeans(randomPresMetrics[,ind]))/apply(randomPresMetrics[,ind],1,sd))
  }, actualPresMetrics, randomPresMetrics)
  Zstatistics = do.call(cbind, Zstatistics)  
  colnames(Zstatistics) = c("meanAdj", "meanClCoeff", "cor.Adj", "cor.kIM", "cor.Cl.Coef")
  Zstatistics = as.data.frame(Zstatistics)
  Zstatistics[,'Z.min'] = apply(Zstatistics[,c(1,3,4)], 1, min)
  Zstatistics[,'Z.median'] = apply(Zstatistics[,c(1,3,4)], 1, median)
  Zstatistics = cbind(actualPresMetrics[,1:2], Zstatistics)
  
  return(list(presMetrics = actualPresMetrics, presRank = actualPresMetrics.rank, randomPresMetrics = randomPresMetrics, Z = Zstatistics))
}
