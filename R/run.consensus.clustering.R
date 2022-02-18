#' Runs Consensus Clustering Algorithm 
#'
#' A modiefed parallel version of code imported from 
#' https://github.com/Bioconductor-mirror/ConsensusClusterPlus 1.11.1. Taken out 
#' of findModules.consensusCluster.R
#'
#' @param d Required. A matrix where columns=items/samples and rows are features. 
#' For example, a gene expression matrix of genes in rows and microarrays in columns.
#' OR ExpressionSet object. (Default = NULL)
#' @param kGrid Optional. (Default = NULL)
#' @param diss Optional. A distance matrix of d. (Default = inherits( d, "dist" ))
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
#' @param weightsItem Optional. A numerical vector. weights to be used for 
#' sampling items. (Default = NULL)
#' @param weightsFeature Optional. AN umerical vector. weights to be used for 
#' sampling features. (Default = NULL)
#' @param verbose Optional. A boolean when set to TRUE, prints messages to the 
#' screen to indicate progress. This is useful for large datasets.(Default = FALSE)
#' @param corUse Optional. Use all cores avaiable. (Default = "Everything")
#'
#' @return Clustered consensus matrix. areaUnderCDF = areaK, 
#' consensus.matrix = cns.mtrx
#' 
#' @importFrom magrittr %>%
#' @export
run.consensus.cluster <- function( d,
                                   kGrid=NULL,
                                   repCount=NULL,
                                   diss=inherits( d, "dist" ),
                                   pItem=NULL,
                                   pFeature=NULL,
                                   innerLinkage=NULL,
                                   distance=NULL, 
                                   clusterAlg=NULL,
                                   weightsItem=NULL,
                                   weightsFeature=NULL,
                                   verbose=NULL,
                                   corUse=NULL) {
  
  n = ifelse( diss, ncol( as.matrix(d) ), ncol(d) )
  
  if (is.null( distance ) ) distance <- 'euclidean'  
  
  acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski", "pearson", "spearman" )
  
  main.dist.obj <- NULL
  if ( diss ){
    main.dist.obj <- d
    ## reset the pFeature & weightsFeature params if they've been set (irrelevant if d is a dist matrix)
    if ( ( !is.null(pFeature) ) && ( pFeature < 1 ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified pFeature parameter\n" )
      pFeature <- 1 # set it to 1 to avoid problems with sampleCols
    }
    if ( ! is.null( weightsFeature ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified weightsFeature parameter\n" )
      weightsFeature <- NULL  # set it to NULL to avoid problems with sampleCols
    }
  } else { ## d is a data matrix
    ## we're not sampling over the features
    if ( ( clusterAlg != "kmeans" ) && ( is.null( pFeature ) || ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) ) {
      # only generate a main.dist.object IFF 
      #    1) d is a distance matrix, 
      #    2) we're not sampling the features, and 
      #    3) the algorithm isn't 'kmeans'
      if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &  ( class(try(get(distance),silent=T))!="function") ) 
          stop("unsupported distance.")
        if(distance=="pearson" | distance=="spearman"){
          main.dist.obj <- stats::as.dist( 1-stats::cor(d,method=distance,use=corUse ))
        } else if( class(try(get(distance),silent=T))=="function"){
          main.dist.obj <- get(distance)( t( d )   )
        } else {
          main.dist.obj <- stats::dist( t(d), method=distance )
        }
        attr( main.dist.obj, "method" ) <- distance  
      } else {
        stop("unsupported distance specified.")
      }
      
    } else {
      ## pFeature < 1 or a weightsFeature != NULL
      ## since d is a data matrix, the user wants to sample over the gene features, so main.dist.obj is left as NULL
    }
  }
  
  cls = plyr::llply(seq(1, repCount, 1), 
                    .fun = function(i, verbose, d, kGrid, pItem, pFeature, weightsItem, weightsFeature, 
                                    main.dist.obj, clusterAlg, distance, acceptable.distance,
                                    corUse, innerLinkage){
                      if(verbose){
                        message(paste("random subsample",i));
                      }
                      
                      # Function to sub sample 
                      sampleCols <- function( d,
                                              pSamp=NULL,
                                              pRow=NULL,
                                              weightsItem=NULL,
                                              weightsFeature=NULL ){
                        ## returns a list with the sample columns, as well as the sub-matrix & sample features (if necessary)
                        ## if no sampling over the features is performed, the submatrix & sample features are returned as NAs
                        ## to reduce memory overhead
                        
                        space <- ifelse( inherits( d, "dist" ), ncol( as.matrix(d) ), ncol(d) )
                        sampleN <- floor(space*pSamp)
                        sampCols <- sort( sample(space, sampleN, replace = FALSE, prob = weightsItem) )
                        
                        this_sample <- sampRows <- NA
                        if ( inherits( d, "matrix" ) ) {
                          if ( (! is.null( pRow ) ) &&
                               ( (pRow < 1 ) || (! is.null( weightsFeature ) ) ) ) {
                            ## only sample the rows and generate a sub-matrix if we're sampling over the row/gene/features
                            space = nrow(d)
                            sampleN = floor(space*pRow)
                            sampRows = sort( sample(space, sampleN, replace = FALSE, prob = weightsFeature) )
                            this_sample <- d[sampRows,sampCols]
                            dimnames(this_sample) <- NULL
                          } else {
                            ## do nothing
                          }
                        }
                        return( list( submat=this_sample,
                                      subrows=sampRows,
                                      subcols=sampCols ) )
                      }
                      
                      # Take expression matrix sample, samples and genes
                      sample_x = sampleCols( d, pItem, pFeature, weightsItem, weightsFeature )
                      
                      # Compute distance (if not supplied)
                      this_dist = NA 
                      if ( ! is.null( main.dist.obj ) ) {
                        boot.cols <- sample_x$subcols
                        this_dist <- as.matrix( main.dist.obj )[ boot.cols, boot.cols ]
                        this_dist <- stats::as.dist( this_dist )
                      } else {
                        # If main.dist.obj is NULL, then d is a data matrix, and either:
                        #   1) clusterAlg is 'kmeans'
                        #   2) pFeatures < 1 or weightsFeatures have been specified, or
                        #   3) both
                        # so we can't use a main distance object and for every iteration, we will have to re-calculate either
                        #   1) the distance matrix (because we're also sampling the features as well), or
                        #   2) the submat (if using kmeans) 
                        
                        if ( clusterAlg != "kmeans" )  {
                          if ( ! distance %in% acceptable.distance &  ( class(try(get(distance),silent=T))!="function")  )
                            stop("unsupported distance.")
                          
                          if( ( class(try(get(distance),silent=T))=="function") ){
                            this_dist <- get(distance)( t( sample_x$submat ) )
                          } else {
                            if( distance == "pearson" | distance == "spearman"){
                              this_dist <- stats::as.dist( 1-stats::cor(sample_x$submat, use=corUse, method=distance) )
                            } else {
                              this_dist <- stats::dist( t( sample_x$submat ), method= distance  )
                            }
                          }
                          attr( this_dist, "method" ) <- distance  
                        } else {
                          # if we're not sampling the features, then grab the colslice
                          if ( is.null( pFeature ) || ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) {
                            this_dist <- d[, sample_x$subcols ]
                          } else {
                            if ( is.na( sample_x$submat ) ) {
                              stop( "error submat is NA" )
                            }
                            this_dist <- sample_x$submat
                          } 
                        }
                      }
                      
                      # Cluster samples using HC (hier. clustering)
                      this_cluster = NA
                      if(clusterAlg == "hc"){
                        this_cluster = fastcluster::hclust(this_dist, method = innerLinkage)
                      }
                      
                      # For every k cluster and identify memebers
                      cls = plyr::llply(kGrid, .fun = function(k, verbose, clusterAlg, this_dist, this_cluster, mConsist, sample_x){
                        if(verbose){
                          message(paste("  k =",k))
                        }
                        
                        this_assignment=NA
                        if(clusterAlg == "hc"){
                          # Prune to k for hc
                          this_assignment = stats::cutree(this_cluster,k)
                        } else if(clusterAlg=="kmeans"){
                          this_assignment = stats::kmeans(t(this_dist), k, iter.max = 100, nstart = 1, algorithm = c("Hartigan-Wong") )$cluster
                        } else if ( clusterAlg == "pam" ) {
                          this_assignment <- cluster::pam(this_dist, k, diss=TRUE, cluster.only=TRUE)
                        } else {
                          # Optional cluterArg Hook.
                          this_assignment <- get(clusterAlg)(this_dist, k)
                        }
                        
                        # Get connectivity matrix				
                        names( this_assignment ) <- sample_x[[3]] 
                        cls <- lapply( unique( this_assignment ), function(clustnum) {
                          as.numeric( names( this_assignment[ this_assignment %in% clustnum ] ) )
                        })  # list samples by clusterId
                        return(cls)
                      },
                      verbose, clusterAlg, this_dist, this_cluster, mConsist, sample_x,
                      .parallel = F,
                      .paropts = list(.packages = c('cluster', 'fastcluster')))
                      names(cls) = kGrid
                      
                      return(cls)
                    },
                    verbose, d, kGrid, pItem, pFeature, weightsItem, weightsFeature, 
                    main.dist.obj, clusterAlg, distance, acceptable.distance,
                    corUse, innerLinkage,
                    .parallel = F,
                    .paropts = list(.packages = c('cluster', 'fastcluster')))
  
  # Compute consensus fraction and area under the cdf curve
  areaK = rep(0,length(cls[[1]]))
  names(areaK) = names(cls[[1]])
  for (k in names(cls[[1]])){
    # Compute consensus matrix
    cns.mtrx = matrix(0, n, n)
    for (nrep in 1:repCount){
      for(nclust in cls[[nrep]][[k]]){
        cns.mtrx[nclust, nclust] = cns.mtrx[nclust, nclust] + 1
      }
    }
    cns.mtrx = cns.mtrx / repCount
    
    # Empirical CDF distribution. default number of breaks is 100    
    h = graphics::hist(cns.mtrx, plot=FALSE, breaks = seq(0,1,by=1/100))
    h$counts = cumsum(h$counts)/sum(h$counts)
    
    # Calculate area under CDF curve, by histogram method.
    thisArea = 0
    for (bi in 1:(length(h$breaks)-1)){
      thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
      bi = bi + 1
    }
    areaK[k] = thisArea
  }
  return(list(areaUnderCDF = areaK, consensus.matrix = cns.mtrx))
}