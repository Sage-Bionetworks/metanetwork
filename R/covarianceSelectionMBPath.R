covarianceSelectionMBPath = function(X,rankedEdges,numberObservations){
  #rankedEdges: list of ranked edges
  nedges <- nrow(rankedEdges)
  count <- 1
  bic <- rep(0,nedges)
  library(utilityFunctions)
  neighborhoods <- vector('list',ncol(X))
  names(neighborhoods)<- colnames(X)
  bicNeighborhood <- apply(X,2,utilityFunctions::fastlmbic)
  names(bicNeighborhood) <- colnames(X)
  bicCurrent <- sum(bicNeighborhood)
  while(count<=nedges){
    print(bicCurrent)
    gene1 <- colnames(X)[rankedEdges[count,1]]
    gene2 <- colnames(X)[rankedEdges[count,2]]
    neighborhoods[[gene1]] <- c(neighborhoods[[gene1]],gene2)
    neighborhoods[[gene2]] <- c(neighborhoods[[gene2]],gene1)
    bicCurrent <- bicCurrent - sum(bicNeighborhood[c(gene1,gene2)])
    bicNeighborhood[gene1] <- utilityFunctions::fastlmbic(X[,gene1],X[,neighborhoods[[gene1]]])
    bicNeighborhood[gene2] <- utilityFunctions::fastlmbic(X[,gene2],X[,neighborhoods[[gene2]]])
    bicCurrent <- bicCurrent + sum(bicNeighborhood[c(gene1,gene2)])
    bic[count] <- bicCurrent
    #plot(bic)
    count <- count+1
  }
  return(bic)
}