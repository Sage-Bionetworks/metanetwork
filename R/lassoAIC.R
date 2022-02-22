#' Grab AIC Solution for Lasso
#' 
#' This function runs a glmnet() function with lasso regression and pulls the best AIC estimate.
#' 
#' @param x Required. An input matrix, of dimension nobs x nvars; each row is an observation vector. 
#' Can be in sparse matrix format (inherit from class "sparseMatrix" as in package Matrix)
#' @param y Required. response variable. Quantitative for family="gaussian", or 
#' family="poisson" (non-negative counts). For family="binomial" should be either 
#' a factor with two levels, or a two-column matrix of counts or proportions 
#' (the second column is treated as the target class; for a factor, the last 
#' level in alphabetical order is the target class). For family="multinomial", 
#' can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. 
#' For either "binomial" or "multinomial", if y is presented as a vector, it will 
#' be coerced into a factor. For family="cox", preferably a Surv object from the 
#' survival package: see Details section for more information. For family="mgaussian",
#'y is a matrix of quantitative responses.
#' 
#' @return Lowest AIC value.
#' 
#' @export
lassoAIC <- function(y,x){
  #require(glmnet)
  res <- glmnet::glmnet(y=y,x=x)
  resid <- y - x%*%res$beta;
  error <- apply(resid^2,2,mean)
  nonzero <- apply(res$beta,2,function(x) sum(x!=0))
  n <- nrow(x)
  bic <- n*log(error)+nonzero*2
  return(res$beta[,which.min(bic)])
}

#bic:nlogsig + klogn