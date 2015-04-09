applyScaleFree <- function(network){
  require(WGCNA)
  network <- network/2 + t(network)/2
  network <- as.matrix(network)
  maxVar <- max(abs(network))
  network <- network/maxVar
  edgeVec <- network[which(upper.tri(network))]
  edgeVec <- abs(edgeVec)
  edgeVec <- edgeVec[which(edgeVec!=0)]
  quant <- (2^(seq(0.5,8,length.out=17))-1)/(2^(seq(0.5,8,length.out=17)))
  cutPoints <- quantile(edgeVec,quant)
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