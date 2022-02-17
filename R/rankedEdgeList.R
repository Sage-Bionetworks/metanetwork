#' Rank edges of a coexpression matrix
#' 
#' This functions produces a ranked edgelist of a coexpression matrix
#' 
#' @param network A network object as a matrix object
#' @param symmetric Optional. A logical `TRUE` or `FALSE` to treat `network` as
#' a symmetrical matrix. (Default = FALSE)
#' @param maxLength Optional. The maximum length of an output ranked edgelist.
#' (Default = 1e7)
#' 
#' @return A greatest to least ranked list of eddge strenge in the form of a 
#' data.frame() object. 
#' 
#' @importFrom magrittr %>%
#' @export
rankedEdgeList <- function(network,symmetric=FALSE,maxLength=1e7){
  #require(dplyr)
  #edgeMat <- matrix(paste0('e',1:(nrow(network)*ncol(network))),nrow(network),ncol(network))
  if(!symmetric){
    #need to fix this
    whichMatrix <- ((network%>%abs) > 0) %>% which(T)
  }else {
    tl <- min(choose(ncol(network),2),maxLength)
    cat('sorting\n')
    a <- sort(abs(network[which(upper.tri(network))]),decreasing=T)[tl]
    a <- max(0,a)
    gc()
    cat('getting matrix\n')
    whichMatrix <- ((network%>%abs) > a & network%>%upper.tri) %>% which(T)    
    gc()
  }
  internal <- function(ind,x){ return(x[ind[1],ind[2]])}
  rankedEdgeList <- cbind(whichMatrix[,1],whichMatrix[,2],apply(whichMatrix,1,internal,network))
  cat('building table\n')
  gc()
  colnames(rankedEdgeList) <- c('var1','var2','value')
  cat('getting edge names\n')
  #rownames(rankedEdgeList) <- apply(whichMatrix,1,internal,edgeMat)
  gc()
  cat('making a dataframe \n')
  rankedEdgeList <- rankedEdgeList %>% data.frame(stringsAsFactors = F)
  gc()
  cat('changing value to numeric\n')
  rankedEdgeList$value <- as.numeric(rankedEdgeList$value)
  gc()
  cat('ordering data frame\n')
  rankedEdgeList <- rankedEdgeList[order((rankedEdgeList$value %>% abs),decreasing=T),]
  gc()
  rankedEdgeList %>% return
}