pickHardThreshold.fromSimilarityMetaNet <- function(similarity, RsquaredCut = 0.85, cutVector = seq(0.1,0.9, by = 0.05), moreNetworkConcepts = FALSE, removeFirst = FALSE,nBreaks = 10) 
{
  checkSimilarity(similarity)
  pickHardThreshold(similarity, dataIsExpr = FALSE, RsquaredCut = RsquaredCut, 
                    cutVector = cutVector, moreNetworkConcepts = moreNetworkConcepts, 
                    removeFirst = removeFirst, nBreaks = nBreaks, corFnc = "I", 
                    corOptions = "")
}