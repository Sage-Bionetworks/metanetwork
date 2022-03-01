#' Applies Bonferroni Correction to a Sparrow Network Object 
#' 
#' This function applies bonferroni correction to a sparrow network object.
#' 
#' @param network A n x n upper triangular adjacency in the matrix class format.
#'
#' @return  A n x n upper triangular adjacency in the matrix class format.
#' 
#' @importFrom magrittr %>%
#' @export
applySparrowBonferroni <- function(network){
  #- require(dplyr) - will add to namespace
  network <- network/2 + t(network)/2
  #network <- as.matrix(network)
  network <- network %>% 
    as.matrix
  #thres <- qchisq(0.05/choose(nrow(network),2),1,lower.tail=F)
  thres <- (0.05 / (nrow(network) %>% choose(2))) %>% stats::qchisq(1, lower.tail=F)
  network <- network^2 > thres
  return(network)
}