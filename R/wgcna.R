# ###Function to run wgcna on the data
library(WGCNA)
wgcna <- function(data,path=NULL){
  softThreshold = pickSoftThreshold(data)$powerEstimate
  network <- abs(cor)^softThreshold
  save(network,file='result_wgcna.rda')
}
