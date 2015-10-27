applySparrowFDR <- function(network){
  require(dplyr)
  network <- network/2 + t(network)/2
  network <- network %>% as.matrix
  #thres <- qchisq(0.05/choose(nrow(network),2),1,lower.tail=F)
  #network1 <- pnorm(abs(network),lower.tail=F)*2
  network1 <- 2*(network %>% abs %>% pnorm(lower.tail=F))
  network1vec <- network1[network1 %>% upper.tri %>% which] %>% c
  thres <- network1vec %>% fdrThres
  network <- network1<thres
  return(network)
}