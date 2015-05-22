rankedEdgeList <- function(network){
  require(dplyr)
  whichMatrix <- ((network%>%abs) > 0) %>% which(T)
  internal <- function(ind,x){ return(x[ind[1],ind[2])}
  cbind(rownames(network)[whichMatrix[,1]],colnames(network)[whichMatrix[,2]],apply(whichMatrix,1,internal,network)) %>% return
}