computeBICcurve <- function(network,exprData,maxEdges=NULL){
  library(dplyr)  
  if(is.null(maxEdges)){
    maxEdges <- round((nrow(exprData)*ncol(exprData))/10)
  }
  
  #networkObj <- synGet('syn5652395',downloadLocation='./')
  #network <- data.table::fread(networkObj@filePath,stringsAsFactors=FALSE,data.table=F)
  #rownames(network) <- network$V1
  #network <- network[,-1]
  #exprData <- read.csv('cranioRNAseq.csv',stringsAsFactors=F,row.names=1)
  
  foo <- data.matrix(network)[which(upper.tri(data.matrix(network)))]
  foo <- foo^2
  thresVal <- sort(foo,decreasing=T)[min(maxEdges,length(foo))]
  
  #add in check for zero edges
  if(thresVal==0){
    thresVal <- min(foo[which(foo>0)])
  }
  
  network <- network^2
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
  
  return(list(bicPath=bicPath,edgeList=edgeList,bicMin=which.min(bicPath$bic)))
}