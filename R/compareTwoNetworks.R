compareTwoNetworks <- function(network1,network2){
  #get networks into the same format
  geneName1 <- colnames(network1)
  geneName2 <- colnames(network2)
  
  geneName <- intersect(geneName1,geneName2)
  
  network1 <- network1[geneName,geneName]
  network2 <- network2[geneName,geneName]
  
  net1 <- network1[which(upper.tri(network1))]
  net2 <- network2[which(upper.tri(netw))]
  
  #compute edge overlaps: OR, p-value from hypergeometric, along with percent overlap
  
  
  
  #compute intersection network
  
  #make degree distribution plots: histogram and scatterplots of degree scores with correlations
    
}