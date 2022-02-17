#' Run genie3
#' 
#' This function runs genie3 random forest network analysis on
#' 
#' @param x Required. A data frame or a matrix of predictors, or a formula 
#' describing the model to be fitted (for the print method, an randomForest object).
#' @param y Required. A response vector. If a factor, classification is assumed,
#' otherwise regression is assumed. If omitted, randomForest will run in 
#' unsupervised mode.
#'
#' @return A matrix of importance measure, one row for each predictor variable. 
#' The column(s) are different importance measures.
#' 
#' @export
genie3 <- function(y,x){
  #function to run genie3 on a single gene
  #IncNodePurity
  #require(randomForest)
  rf <- randomForest::randomForest(x,
                                   y,
                                   mtry=round(sqrt(ncol(x))),
                                   ntree=1000,
                                   importance=TRUE
                                  )
  im <- randomForest::importance(rf)[,'IncNodePurity'];
  return(im)
}
