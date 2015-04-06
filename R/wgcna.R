# ###Function to run wgcna on the data
wgcna <- function(data,path=NULL){
  library(WGCNA)
  softThresholdMatrix = pickSoftThreshold(data)
  if(!is.na(softThresholdMatrix$powerEstimate)){
    softThreshold = softThresholdMatrix$powerEstimate
  }else{
    w1 <- which.max(softThresholdMatrix$fitIndices$SFT.R.sq)
    softThreshold = softThresholdMatrix$fitIndices$Power[w1]
    cat('warning, r^2 less than .85:',softThresholdMatrix$fitIndices$SFT.R.sq[w1],'\n')
  }
  network <- abs(cor(data))^softThreshold
  save(network,file='result_wgcna.rda')
}
