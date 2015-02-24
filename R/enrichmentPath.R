enrichmentPath <- function(targetList,rankedList){
  doEnrich<-function(i,targetList,rankedList){
    foo <- unlist(enrichment(targetList,rankedList[1:i],rankedList))
    foo <- c(rankedList[i],foo)
    names(foo) <- c('geneId','enr','pval')
    return(foo)
  }
  return(data.frame(t(sapply(1:length(rankedList),doEnrich,targetList,rankedList))))
}