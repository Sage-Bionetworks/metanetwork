computeDriverDistance  = function(geneSet,graph){
  library(igraph)
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
  #bar[!is.finite(bar)]<-NA
  return(list(meanScore=apply(bar,1,mean),distMat=bar))
}