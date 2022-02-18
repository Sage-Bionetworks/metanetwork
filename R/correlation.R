#' Compute a Correlation Matrix
#' 
#' Compute a Correlation matrix from a an expression matrix and writes to a
#' user specified path. Correlation network is saved as correlationNetwork.csv on
#' the path specified by `outputpath`.
#' 
#' @param data The input expression/numerical data matrix
#' @param path Optional. Currently not implemented in the function
#' @param method Optional. The type of correlation to be implemented. Options are
#' one of; c('pearson','spearman','kendall'). (Default = pearson)
#' @param outputpath The file path to save the correlation network to.
#' 
#' @return Returns no value direct to R. This function saves the correlation network
#' as correlationNetwork.csv on the path specified by `outputpath`.
#' 
#' @export
#' 
correlation <- function(data,path=NULL,method='pearson',outputpath){
  network <- stats::cor(data,method=method)
  #save(network,file=paste0(outputpath,'result_correlation.rda'))
  network <- network*upper.tri(network)
  utils::write.csv(network,file=paste0(outputpath,'correlationNetwork.csv'),quote=F)
}