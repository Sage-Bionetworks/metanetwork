applyScaleFree <- function(network){
  require(WGCNA)
  require(dplyr)
  network <- network %>% as.matrix
  #maxVar <- max(abs(network))
  maxVar <- (network %>% abs) %>% max
  network <- network/maxVar
  #edgeVec <- network[which(upper.tri(network))]
  edgeVec <- network[(network %>% upper.tri) %>% which] %>% abs
  #edgeVec <- abs(edgeVec)
  edgeVec <- edgeVec[(edgeVec!=0) %>% which]
  quant <- (2^(seq(5,15,length.out=17))-1)/(2^(seq(5,15,length.out=17)))
  #cutPoints <- quantile(edgeVec,quant)
  cutPoints <- edgeVec %>% quantile(quant)
  diag(network) <- 1
  #network <- abs(network)
  network <- network %>% abs
  hardThresholdMatrix <- pickHardThreshold.fromSimilarityMetaNet(similarity = network, cutVector = cutPoints)
  if(!is.na(hardThresholdMatrix$cutEstimate)){
    hardThreshold = hardThresholdMatrix$cutEstimate
  }else{
    w1 <- which.max(hardThresholdMatrix$fitIndices$SFT.R.sq)
    hardThreshold = hardThresholdMatrix$fitIndices$Cut[w1]
    cat('warning, r^2 less than .85:',hardThresholdMatrix$fitIndices$SFT.R.sq[w1],'\n')
  }
  #diag(network) <- 0;
  
  #network <- abs(network) > hardThreshold
  n1<-((network[(network %>% upper.tri) %>% which] %>% abs) > hardThreshold) %>% sum
  return(n1)
}