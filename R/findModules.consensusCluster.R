#' Finds Consensus Clusters
#'
#' A modiefed parallel version of code imported from 
#' https://github.com/Bioconductor-mirror/ConsensusClusterPlus 1.11.1
#'
#' @param d Optional. A matrix where columns=items/samples and rows are features. 
#' For example, a gene expression matrix of genes in rows and microarrays in columns.
#' OR ExpressionSet object. (Default = NULL)
#' @param maxK Optional. An  integer value. maximum cluster number to evaluate.
#' (Default = 100)
#' @param reps Optional. An integer value. number of subsamples.  (Default = 100)
#' @param pItem Optional. A numerical value. proportion of items to sample.
#'  (Default = 0.8)
#' @param pFeature Optional. A numerical value. proportion of features to sample.
#' (Default = 1)
#' @param clusterAlg Optional. A character value. cluster algorithm. "hc" 
#' heirarchical (hclust) or "km" for kmeans. (Default = "kmeans") 
#' @param innerLinkage Optional. A heirarchical linkage method for subsampling.
#' (Default = "average")
#' @param distance Optional. A character value. sample distance measures: 
#' "pearson","spearman", or "euclidean". (Default = "pearson")
#' @param seed Optional, A numerical value. Sets random seed for reproducible results.
#' (Default = 123456789.12345)
#' @param weightsItem Optional. A numerical vector. weights to be used for 
#' sampling items. (Default = NULL)
#' @param weightsFeature Optional. AN umerical vector. weights to be used for 
#' sampling features. (Default = NULL)
#' @param verbose Optional. A boolean when set to TRUE, prints messages to the 
#' screen to indicate progress. This is useful for large datasets.(Default = FALSE)
#' @changeCDFArea Optional. Minimum spline distance for seq(2,`maxK`,
#' length.out = `nbreaks`) (Default = 0.001)
#' @param corUse Optional. Use all cores avaiable. (Default = "Everything")
#' @param nbreaks Optional. Number of breaks to use in 
#' seq(2,`maxK`,length.out = `nbreaks`) this becomes the kGrid argument in 
#' run.consensus.cluster. (Default = 20)
#' 
#' @return  Final clustered modules.
#' 
#' @importFrom magrittr %>%
#' @export
findModules.consensusCluster <- function(d = NULL,
                                         maxK = 100,
                                         reps = 100,
                                         pItem = 0.8,
                                         pFeature = 1,
                                         clusterAlg = "kmeans",
                                         innerLinkage = "average",
                                         distance = "pearson",
                                         changeCDFArea = 0.001,
                                         nbreaks = 20,
                                         seed = 123456789.12345,
                                         weightsItem = NULL,
                                         weightsFeature = NULL,
                                         corUse = "everything",
                                         verbose = F) {
  # Set seed 
  if(is.null(seed)==TRUE){
    seed=timeSeed = as.numeric(Sys.time())
  }
  set.seed(seed)
  
  # Error checking
  if ( ! class( d ) %in% c( "dist", "matrix", "ExpressionSet" ) ) {
    stop("d must be a matrix, distance object or ExpressionSet (eset object)")
  }
    
  # If distance matrix is supplied instead of adj
  if ( inherits( d, "dist" ) ) {
    ## if d is a distance matrix, fix a few things so that they don't cause problems with the analysis
    ##  Note, assumption is that if d is a distance matrix, the user doesn't want to sample over the row features
    if ( is.null( attr( d, "method" ) ) ) {
      attr( d, "method" ) <- distance <- "unknown - user-specified"
    }
    
    if ( is.null( distance ) || ( distance != attr( d, "method" ) ) ) {
      distance <- attr( d, "method" ) 
    }
      
    if ( ( ! is.null( pFeature ) ) && ( pFeature < 1 ) ) {
      message( "Cannot use the pFeatures parameter when specifying a distance matrix as the data object\n" )
      pFeature <- 1
    }
    
    if ( ! is.null( weightsFeature ) ) {
      message( "Cannot use the weightsFeature parameter when specifying a distance matrix as the data object\n" )
      weightsFeature <- NULL
    }
    
    if ( clusterAlg == "kmeans" )
      message( "Note: k-means will cluster the distance matrix you provided.  This is similar to kmdist option when supplying a data matrix")
    
  } else {
    if ( is.null( distance ) ) {
      ## we should never get here, but just in case
      distance <- "pearson"
    }
  }
  
  # If d is a expression set
  if ( inherits( d,"ExpressionSet" ) )
    d <- exprs(d)
  
  # Run consensus clustering
  kGrid = seq(2,maxK,length.out = nbreaks) %>%
    round() %>% unique()
  results <- run.consensus.cluster( d=d,
                                    kGrid=kGrid,
                                    repCount=reps,
                                    diss=inherits(d,"dist"),
                                    pItem=pItem,
                                    pFeature=pFeature,
                                    innerLinkage=innerLinkage,
                                    clusterAlg=clusterAlg,
                                    weightsFeature=weightsFeature,
                                    weightsItem=weightsItem,
                                    distance=distance,
                                    verbose=verbose,
                                    corUse=corUse )
    
  # Check if ends are maximum 
  fn = splinefun(names(results$areaUnderCDF), results$areaUnderCDF)
  ind = which(diff(fn(2:maxK)) >= changeCDFArea)
  k.final = ind[length(ind)] + 1
  
  # Re-run the consensus for k.final
  results.final <- run.consensus.cluster( d=d,
                                          kGrid=k.final,
                                          repCount=reps,
                                          diss=inherits(d,"dist"),
                                          pItem=pItem,
                                          pFeature=pFeature,
                                          innerLinkage=innerLinkage,
                                          clusterAlg=clusterAlg,
                                          weightsFeature=weightsFeature,
                                          weightsItem=weightsItem,
                                          distance=distance,
                                          verbose=verbose,
                                          corUse=corUse )
  
  # Compute the final clusters
  cluster.final = kmeans(1 - results.final$consensus.matrix, k.final, algorithm = "Hartigan-Wong", trace = TRUE)$cluster
  mod = data.frame(Gene.ID = colnames(d),
                   moduleNumber = cluster.final)
  
  # Change cluster number to color labels
  mod$moduleLabel = WGCNA::labels2colors(mod$moduleNumber)
  
  return(mod)
}
