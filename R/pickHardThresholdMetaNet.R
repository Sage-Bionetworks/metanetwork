pickHardThresholdMetaNet <- function(data, dataIsExpr = TRUE, RsquaredCut = 0.85, cutVector = seq(0.1, 
                                                                       0.9, by = 0.05), moreNetworkConcepts = FALSE, removeFirst = FALSE, 
          nBreaks = 10, corFnc = "cor", corOptions = "use = 'p'") 
{
  nGenes = dim(data)[[2]]
  colname1 = c("Cut", "p-value", "SFT.R.sq", "slope=", "truncated R^2", 
               "mean(k)", "median(k)", "max(k)")
  if (moreNetworkConcepts) {
    colname1 = c(colname1, "Density", "Centralization", "Heterogeneity")
  }
  if (!dataIsExpr) {
    checkAdjMat(data)
    if (any(diag(data) != 1)) 
      diag(data) = 1
  }
  else nSamples = dim(data)[[1]]
  datout = data.frame(matrix(NA, nrow = length(cutVector), 
                             ncol = length(colname1)))
  names(datout) = colname1
  datout[, 1] = cutVector
  for (i in 1:length(cutVector)) {
    cut1 = cutVector[i]
    datout[i, 2] = 2 * (1 - pt(sqrt(nSamples - 1) * cut1/sqrt(1 - 
                                                                cut1^2), nSamples - 1))
  }
  if (exists("fun1")) 
    rm(fun1)
  fun1 = function(x, dataIsExpr) {
    if (dataIsExpr) {
      corExpr = parse(text = paste(corFnc, "(x, data", 
                                   prepComma(corOptions), ")"))
      corx = abs(eval(corExpr))
    }
    else corx = x
    out1 = rep(NA, length(cutVector))
    for (j in c(1:length(cutVector))) {
      out1[j] = sum(corx > cutVector[j])
    }
    out1
  }
  datk = t(apply(data, 2, fun1, dataIsExpr))
  for (i in c(1:length(cutVector))) {
    khelp = datk[, i] - 1
    SFT1 = scaleFreeFitIndex(k = khelp, nBreaks = nBreaks, 
                             removeFirst = removeFirst)
    datout[i, 3] = SFT1$Rsquared.SFT
    datout[i, 4] = SFT1$slope.SFT
    datout[i, 5] = SFT1$truncatedExponentialAdjRsquared
    datout[i, 6] = mean(khelp, na.rm = T)
    datout[i, 7] = median(khelp, na.rm = T)
    datout[i, 8] = max(khelp, na.rm = T)
    if (moreNetworkConcepts) {
      Density = sum(khelp)/(nGenes * (nGenes - 1))
      datout[i, 9] = Density
      Centralization = nGenes * (max(khelp) - mean(khelp))/((nGenes - 
                                                               1) * (nGenes - 2))
      datout[i, 10] = Centralization
      Heterogeneity = sqrt(nGenes * sum(khelp^2)/sum(khelp)^2 - 
                             1)
      datout[i, 11] = Heterogeneity
    }
  }
  print(signif(data.frame(datout), 3))
  ind1 = datout[, 3] > RsquaredCut
  indcut = NA
  indcut = ifelse(sum(ind1) > 0, min(c(1:length(ind1))[ind1]), 
                  indcut)
  cutEstimate = cutVector[indcut][[1]]
  list(cutEstimate = cutEstimate, fitIndices = data.frame(datout))
}