overlapStats <- function (labels1, labels2, na.rm = TRUE, ignore = NULL, levels1 = NULL, levels2 = NULL) {
  
  labels1 = as.vector(labels1)
  labels2 = as.vector(labels2)
  
  if (na.rm) {
    keep = !is.na(labels1) & !is.na(labels2)
    labels1 = labels1[keep]
    labels2 = labels2[keep]
  }
  
  if (is.null(levels1)) {
    levels1 = sort(unique(labels1))
    levels1 = levels1[!levels1 %in% ignore]
  }
  
  if (is.null(levels2)) {
    levels2 = sort(unique(labels2))
    levels2 = levels2[!levels2 %in% ignore]
  }
  
  n1 = length(levels1)
  n2 = length(levels2)
  
  countMat = matrix(0, n1, n2)
  pMat = matrix(0, n1, n2)
  oddsRatio = matrix(0, n1, n2)
  
  for (m1 in 1:n1) for (m2 in 1:n2) {
    m1Members = (labels1 == levels1[m1])
    m2Members = (labels2 == levels2[m2])
    tmp = fisher.test(m1Members, m2Members, alternative = "greater")
    pMat[m1, m2] = tmp$p.value
    oddsRatio[m1, m2] = tmp$estimate
    countMat[m1, m2] = sum(labels1 == levels1[m1] & labels2 == levels2[m2])    
  }
  dimnames(pMat) = list(levels1, levels2)
  dimnames(countMat) = list(levels1, levels2)
  dimnames(oddsRatio) = list(levels1, levels2)
  
  pMat[is.na(pMat)] = 1
  
  list(countTable = countMat, pTable = pMat, orTable = oddsRatio)
}
