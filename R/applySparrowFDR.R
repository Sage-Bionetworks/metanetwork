applySparrowFDR <- function(network){
  network <- network/2 + t(network)/2
  network <- as.matrix(network)
  #thres <- qchisq(0.05/choose(nrow(network),2),1,lower.tail=F)
  #network1 <- pnorm(abs(network),lower.tail=F)*2
  require(dplyr)
  network1 <- 2*(network %>% abs() %>% pnorm(lower.tail=F))
  network1vec <- c(network1[which(upper.tri(network1))])
  thres <- fdrThres(network1vec)
  network <- network1<thres
  return(network)
}