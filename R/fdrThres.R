fdrThres <- function (pval, fdr = 0.05) {
  require(dplyr)
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