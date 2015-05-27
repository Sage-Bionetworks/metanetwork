correlation <- function(data,path=NULL,method='pearson',outputpath){
  network <- cor(data,method='pearson')
  save(network,file=paste0(outputpath,'result_correlation.rda'))
}