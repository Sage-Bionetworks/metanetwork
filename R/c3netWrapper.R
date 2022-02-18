#' This function wraps c3net implementation
#' 
#' This function wraps the c3net call on a gene expression matrix. This takes the
#' transpose of typical gene expression matrix where rows are sample IDs and 
#' columns are Gene IDs in order to convert it for calling the `c3net()` function.
#' 
#' @param data Required. A Matrix containing gene expression values with subjects 
#' as row values and gene features as column IDs.
#' @param outputpath Required. Path to save the resulting network
#' @param pval Optional. Currently not implemented. (Default = 1)
#' 
#' @export 
#' @return Returns a symmetric mutual information matrix, which is obtained after
#'  implementing C3NET. Specifically, non-zero elements in the returned matrix 
#' represents undirected link between variables. The inferred network may also 
#' be plotted if the argument network is set TRUE.
#' 
#' @references G. Altay, F. Emmert-Streib, "Inferring the conservative causal core of gene regulatory networks", BMC Systems Biology (2010) 4:132.
#'
c3netWrapper  = function(data,pval=1,outputpath){
  #library(c3net)
  network <- c3net::c3net(t(data))
  #save(network,file=paste0(outputpath,'result_c3net.rda'))
  network <- network*upper.tri(network)
  utils::write.csv(network,file=paste0(outputpath,'c3netNetwork.csv'),quote=F)
}