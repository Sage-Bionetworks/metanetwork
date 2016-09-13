computeDriverDistancePvalue  = function(geneSet,graph,nsamp=100){
  library(igraph)
  library(metanetwork)
  #####extract the largest connected component
  compLabel <- igraph::components(graph)
  maxComp <- which.max(compLabel$csize)
  cat('largest component size: ',max(compLabel$csize),'\n')
  graph <- subgraph(graph, which(compLabel$membership==maxComp))
  
  geneSetIdx <- which(geneSet%in%names(V(graph)))
  if(length(geneSetIdx)<4){
    stop('target set is too small (<4)')
  }else{
    geneSet <- geneSet[geneSetIdx]
  }
  
  bar <- igraph::distances(graph,names(V(graph)),geneSet)
  meanScore <- apply(bar,1,mean)
  nullVec <- rep(0,nsamp)
  for(i in 1:nsamp){
    nullVec[i] <- igraph::distances(graph,names(V(graph)),sample(names(V(graph)),length(geneSet))) %>% apply(1,mean) %>% min
    cat(min(nullVec[1:i]),'\n')
  }
  return(list(distMat=bar,meanScore=meanScore,nullDist=nullVec))
  
  
  
  #bar[!is.finite(bar)]<-NA
  #return(list(meanScore=apply(bar,1,mean),distMat=bar))
}