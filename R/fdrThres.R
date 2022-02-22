#' FDR Threshold
#' 
#' This function applies a user FDR threshold to input p-values
#' 
#' @param pval Required. A vector of uncorected P-Values.
#' @param fdr Optional. desired FDR cutoff. (Default = 0.05)
#' as y
#' @return Corrected PValues
#' 
#' @importFrom magrittr %>%
#' @export
fdrThres <- function (pval, fdr = 0.05) {
  #require(dplyr)
  n <- pval %>% length
  comp <- (pval %>% sort) < ((fdr/n) * (1:n))
  if (min(pval) < (fdr/n)) {
    w1 <- which(!comp)[1]
    return((fdr/n) * w1)
  }
  else {
    return(fdr/n)
  }
}