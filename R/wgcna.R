# ###Function to run wgcna on the data
wgcna <- function(data,path=NULL){
  library(WGCNA)
  softThreshold = pickSoftThreshold(data)$powerEstimate
  network <- abs(cor)^softThreshold
  save(network,file='result_wgcna.rda')
}
