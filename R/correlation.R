correlation <- function(data,path=NULL,method='pearson'){
  network <- cor(data,method='pearson')
  save(network,file='result_correlation.rda')
}