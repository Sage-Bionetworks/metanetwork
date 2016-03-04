computeBICcurve <- function(network,exprData,maxEdges=NULL,exact=NULL){
  library(dplyr)
  if(is.null(maxEdges)){
    maxEdges <- round((nrow(exprData)*ncol(exprData))/25)
  }

  foo <- data.matrix(network)[which(upper.tri(data.matrix(network)))]
  foo <- abs(foo)
  network[which(lower.tri)]<-0
  diag(network) <- 0
  thresVal <- sort(foo,decreasing=T)[min(maxEdges,length(foo))]

  #add in check for zero edges
  if(thresVal==0){
    thresVal <- min(foo[which(foo>0)])
  }

  network <- abs(network)
  edgeList <- which(network >=thresVal,T)
  edval <- rep(0,nrow(edgeList))

  for (i in 1:nrow(edgeList)){
    edval[i] <- network[edgeList[i,1],edgeList[i,2]]
  }

  edgeList <- cbind(edgeList,edval)
  colnames(edgeList) <- c('node1','node2','weight')
  rownames(edgeList) <- paste0('e',1:nrow(edgeList))
  edgeList <- data.frame(edgeList,stringsAsFactors=F)
  edgeList <- arrange(edgeList,desc(weight))
  bicPath <- metanetwork::covarianceSelectionMBPath(data.matrix(exprData),rankedEdges=edgeList[,1:2],start=1)
  bicPath2 <- NA;

  #if(!is.null(exact)){
  #  bicPath2 <- metanetwork::covarianceSelectionPath(cor(exprData),edgeList[,1:2],nrow(exprData)*ncol(exprData),nrow(exprData))
  #}
  library(Matrix)
  network <- network>=edgeList$weight[which.min(bicPath$bic)]
  return(list(network=Matrix(network,sparse=T),bicMin=which.min(bicPath$bic)))
}
