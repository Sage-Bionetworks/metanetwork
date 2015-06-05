rankedEdgeList <- function(network,symmetric=FALSE,maxLength=1e7){
  require(dplyr)
  edgeMat <- matrix(paste0('e',1:(nrow(network)*ncol(network))),nrow(network),ncol(network))
  #rownames(edgeMat) <- rownames(network)
  #colnames(edgeMat) <- colnames(netw)
  if(!symmetric){
    #need to fix this
    whichMatrix <- ((network%>%abs) > 0) %>% which(T)
  }else {
    tl <- min(choose(ncol(network),2),maxLength)
    a <- sort(abs(network[which(upper.tri(network))]),decreasing=T)[tl]
    a <- max(0,a)
    gc()
    #print(a)
    whichMatrix <- ((network%>%abs) > a & network%>%upper.tri) %>% which(T)    
  }
  internal <- function(ind,x){ return(x[ind[1],ind[2]])}
  rankedEdgeList <- cbind(whichMatrix[,1],whichMatrix[,2],apply(whichMatrix,1,internal,network))
  gc()
  colnames(rankedEdgeList) <- c('var1','var2','value')
  #rownames(rankedEdgeList) <- paste0('edge',1:nrow(rankedEdgeList))
  rownames(rankedEdgeList) <- apply(whichMatrix,1,internal,edgeMat)
  gc()
  rankedEdgeList <- rankedEdgeList %>% data.frame(stringsAsFactors = F)
  gc()
  rankedEdgeList$value <- as.numeric(rankedEdgeList$value)
  gc()
  rankedEdgeList <- rankedEdgeList[order((rankedEdgeList$value %>% abs),decreasing=T),]
  gc()
  #rownames(rankedEdgeList) <- paste0('edge',1:nrow(rankedEdgeList))
  rankedEdgeList %>% return
}