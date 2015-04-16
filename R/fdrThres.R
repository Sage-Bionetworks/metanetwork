fdrThres <- function (pval, fdr = 0.05) {
  n <- length(pval)
  comp <- sort(pval) < ((fdr/n) * (1:n))
  if (min(pval) < (fdr/n)) {
    w1 <- which(!comp)[1]
    return((fdr/n) * w1)
  }
  else {
    return(fdr/n)
  }
}