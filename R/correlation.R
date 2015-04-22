correlation <- function(data,method='pearson'){
  network <- cor(data,method='pearson')
  save(network,file='result_correlation.rda')
}