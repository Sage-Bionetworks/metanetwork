#' Applies FDR Correction to a Sparrow Network Object 
#' 
#' This function applies FDR correction to a sparrow network object.
#' 
#' @param network A n x n upper triangular adjacency in the matrix class format.
#'
#' @return  A n x n upper triangular adjacency in the matrix class format.
#' 
#' @importFrom magrittr %>%
#' @export
applySparrowFDR <- function(network){
  # require(dplyr)
  network <- network/2 + t(network)/2
  network <- network %>% as.matrix
  #thres <- qchisq(0.05/choose(nrow(network),2),1,lower.tail=F)
  #network1 <- pnorm(abs(network),lower.tail=F)*2
  network1 <- 2*(network %>% abs %>% stats::pnorm(lower.tail=F))
  network1vec <- network1[network1 %>% upper.tri %>% which] %>% c
  thres <- network1vec %>% fdrThres
  network <- network1<thres
  return(network)
}