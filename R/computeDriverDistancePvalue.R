#' Compute Driver Distance P-Value
#'
#' Calculates the Driver Distance of `geneSet` in `graph` as described in
#' `metanetworks::computeDriverDistance` and permutes `nsamp` permutations of random
#' gene sets from `graph` which are the size of `geneSet` to create a null distribution
#' to derive p-values for the Driver Distances of `geneSet`.
#'
#' @param nsamp Optional. Number of permutations to build null distribution for
#' p-value construction (Default = 100)
#' 
#' @inheritParams computeDriverDistance
#' 
#' @return a list consisting of a mean distance of genes within the largest component 
#' to each other, the pairwise list of distances, and a vector the length of `nsamp`
#' consisting of the null distribution to calculate the p-value.
#'
#' @importFrom magrittr %>%
#' @export
computeDriverDistancePvalue  = function(geneSet,graph,nsamp=100){
  #library(igraph)
  #library(metanetwork)
  #####extract the largest connected component
  compLabel <- igraph::components(graph)
  maxComp <- which.max(compLabel$csize)
  cat('largest component size: ',max(compLabel$csize),'\n')
  graph <- igraph::subgraph(graph, which(compLabel$membership==maxComp))
  
  geneSetIdx <- which(geneSet%in%names(igraph::V(graph)))
  if(length(geneSetIdx)<4){
    stop('target set is too small (<4)')
  }else{
    geneSet <- geneSet[geneSetIdx]
  }
  
  bar <- igraph::distances(graph,names(igraph::V(graph)),geneSet)
  meanScore <- apply(bar,1,mean)
  nullVec <- rep(0,nsamp)
  for(i in 1:nsamp){
    nullVec[i] <- igraph::distances(graph,
                                    names(igraph::V(graph)),
                                    sample(
                                      names(igraph::V(graph)),
                                      length(geneSet))
                                    ) %>% 
      apply(1,mean) %>% 
      min
    
    cat(min(nullVec[1:i]),'\n')
  }
  return(list(distMat=bar,meanScore=meanScore,nullDist=nullVec))
  #bar[!is.finite(bar)]<-NA
  #return(list(meanScore=apply(bar,1,mean),distMat=bar))
}