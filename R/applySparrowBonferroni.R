applySparrowBonferroni <- function(network){
  require(dplyr)
  network <- network/2 + t(network)/2
  #network <- as.matrix(network)
  network <- network %>% as.matrix
  #thres <- qchisq(0.05/choose(nrow(network),2),1,lower.tail=F)
  thres <- (0.05/(nrow(network) %>% choose(2))) %>% qchisq(1,lower.tail=F)
  network <- network^2>thres
  return(network)
}