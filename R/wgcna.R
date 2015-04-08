# ###Function to run wgcna on the data
wgcna <- function(data,path=NULL){
  library(WGCNA)
  hardThresholdMatrix = pickHardThreshold(data)
  if(!is.na(hardThresholdMatrix$cutEstimate)){
    hardThreshold = hardThresholdMatrix$cutEstimate
  }else{
    w1 <- which.max(hardThresholdMatrix$fitIndices$SFT.R.sq)
    hardThreshold = hardThresholdMatrix$fitIndices$Cut[w1]
    cat('warning, r^2 less than .85:',hardThresholdMatrix$fitIndices$SFT.R.sq[w1],'\n')
  }
  network <- abs(cor(data))>hardThreshold
  save(network,file='result_wgcna.rda')
}
