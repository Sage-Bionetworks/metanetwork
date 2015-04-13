makeSparseNetworkFile <- function(network,geneId,file){
  network <- network*upper.tri(network)
  w1 <- which(network==1,arr.ind = TRUE)
  mat <- cbind(geneId[w1[,1]],geneId[w1[,2]])
  write.csv(mat,file=file,row.names=F,col.names=F)
}