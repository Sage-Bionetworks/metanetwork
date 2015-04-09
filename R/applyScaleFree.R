applyScaleFree <- function(network){
  require(WGCNA)
  network <- network/2 + t(network)/2
  network <- as.matrix(network)
  maxVar <- max(abs(network))
  network <- network/maxVar
  edgeVec <- network[which(upper.tri(network))]
  edgeVec <- abs(edgeVec)
  edgeVec <- edgeVec[which(edgeVec!=0)]
  cutPoints <- quantile(edgeVec,seq(0.91,0.99,length.out=17))
  diag(network) <- 1
  hardThresholdMatrix <- pickHardThreshold.fromSimilarityMetaNet(similarity = network, cutVector = cutPoints)
  if(!is.na(hardThresholdMatrix$cutEstimate)){
    hardThreshold = hardThresholdMatrix$cutEstimate
  }else{
    w1 <- which.max(hardThresholdMatrix$fitIndices$SFT.R.sq)
    hardThreshold = hardThresholdMatrix$fitIndices$Cut[w1]
    cat('warning, r^2 less than .85:',hardThresholdMatrix$fitIndices$SFT.R.sq[w1],'\n')
  }
  network <- abs(network) > hardThreshold
  return(network)
}