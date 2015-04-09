applySparrowBonferroni <- function(network){
  network <- network/2 + t(network)/2
  network <- as.matrix(network)
  thres <- qchisq(0.05/choose(nrow(network),2),1,lower.tail=F)
  network <- network^2>thres
  return(network)
}