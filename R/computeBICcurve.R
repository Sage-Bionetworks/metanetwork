computeBICcurve <- function(network,exprData,maxEdges=NULL,exact=NULL){
  library(dplyr)
  if(is.null(maxEdges)){
    maxEdges <- round((nrow(exprData)*ncol(exprData))/20)
  }
  maxEdges <- min(maxEdges,round((nrow(exprData)*ncol(exprData))/20)) 
  cat('maxEdges:',maxEdges,'\n')
  foo <- data.matrix(network)[which(upper.tri(data.matrix(network)))]
  cat('organize foo\n')
  foo <- abs(foo)
  cat('set lower triangular to 0\n')
  network[(lower.tri(network))==TRUE]<-0
  cat('set diagonal to zero\n')
  diag(network) <- 0
 
  thresVal <- sort(foo,decreasing=T)[min(maxEdges,length(foo))]
  cat('thresVal:',thresVal,'\n')
  #add in check for zero edges
  if(thresVal==0){
    thresVal <- min(foo[which(foo>0)])
  }
  cat('network absolution\n')
  network <- abs(network)
  cat('get edges\n')
  edgeList <- which(network >=thresVal,T)
  edval <- rep(0,nrow(edgeList))

  for (i in 1:nrow(edgeList)){
    if(i%%1000==0){
      cat(edval[i],'\n')
    }
    edval[i] <- network[edgeList[i,1],edgeList[i,2]]
  }

  edgeList <- cbind(edgeList,edval)
  colnames(edgeList) <- c('node1','node2','weight')
  rownames(edgeList) <- paste0('e',1:nrow(edgeList))
  edgeList <- data.frame(edgeList,stringsAsFactors=F)
  edgeList <- arrange(edgeList,desc(weight))
  print(edgeList[1:5,])
  bicPath <- metanetwork::covarianceSelectionMBPath(data.matrix(exprData),rankedEdges=edgeList[,1:2],start=1)
  bicPath2 <- NA;

  #if(!is.null(exact)){
  #  bicPath2 <- metanetwork::covarianceSelectionPath(cor(exprData),edgeList[,1:2],nrow(exprData)*ncol(exprData),nrow(exprData))
  #}
  library(Matrix)
  network <- network>=edgeList$weight[which.min(bicPath$bic)]
  return(list(network=Matrix(network,sparse=T),bicMin=min(bicPath$bic),bicPath=bicPath$bic))
}
