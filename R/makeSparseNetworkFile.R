makeSparseNetworkFile <- function(network,geneId,file){
  network <- network*upper.tri(network)
  w1 <- which(network==1,arr.ind = TRUE)
  mat <- cbind(geneId[w1[,1]],geneId[w1[,2]])
  write.table(mat,file=file,row.names=FALSE,col.names=FALSE,quote=FALSE,sep=',')
}