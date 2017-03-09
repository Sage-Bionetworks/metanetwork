compareTwoNetworks <- function(network1,network2){
  #get networks into the same format
  geneName1 <- colnames(network1)
  geneName2 <- colnames(network2)
  geneName <- intersect(geneName1,geneName2)
  network1 <- network1[geneName,geneName]
  network2 <- network2[geneName,geneName]
  net1 <- network1[which(upper.tri(network1))]
  net2 <- network2[which(upper.tri(network2))]
  table1 <- table(net1!=0,net2!=0)
  model <- list()
  model$fisher <- fisher.test(table1)
  model$overlapNet1 <- table1[1,1]/(rowSums(table1)[1])
  model$overlapNet2 <- table1[2,2]/(rowSums(table1)[2])
  return(model)
}