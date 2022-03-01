#' Compute Graph Distance
#'
#' Finds the finds the maximal (weakly or strongly) connected component of `graph` 
#' then returns the distances between members in `geneSet` which are found in the 
#' largest component.
#'
#' @param geneSet Required. A user specified gene set corresponding to vertex names
#' in `graph`. All the names are not required to be entirely represented in the 
#' graph object.
#' @param graph Required. An igraph graph object consisting of vertices of genes
#' and edges representing co-expression between genes.
#'
#' @return a list consisting of a mean distance of genes within the largest component 
#' to each other and the paiwise list of distances.
#'
#' @export
computeDriverDistance  = function(geneSet,graph){
  #library(igraph)
  #####extract the largest connected component
  compLabel <- igraph::components(graph)
  maxComp <- which.max(compLabel$csize)
  cat('largest component size: ',max(compLabel$csize),'\n')
  graph <- igraph::subgraph(graph, which(compLabel$membership==maxComp))
  
  geneSetIdx <- which(geneSet %in% names(igraph::V(graph)))
  if(length(geneSetIdx)<4){
    stop('target set is too small (<4)')
  }else{
    geneSet <- geneSet[geneSetIdx]
  }
  
  bar <- igraph::distances(graph, names(igraph::V(graph)), geneSet)
  #bar[!is.finite(bar)]<-NA
  return(list(meanScore=apply(bar,1,mean),distMat=bar))
}