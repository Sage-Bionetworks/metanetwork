correlation <- function(data){
  network <- cor(data)
  save(network,file='result_correlation.rda')
}