#' Runs vbsr gene across a matrix 
#' 
#' This function wraps variable bays spike regression of a genes expression across
#' a matrix of genes expressed in the same samples.
#' 
#' @param y Required. response variable. Normally distributed errors for 
#' family="normal". For family="binomial" should be coded as a vector of 0's and 1's.
#' @param x Required. Design matrix, an n x m matrix, with rows as observations.
#' @param ordering_mat Optional. Optionally specified coordinate update ordering matrix.
#'  Must be in matrix form with columns as permutation vectors of length m, and 
#'  there must be n_orderings columns. (Default = NULL)
#' @param eps Optional. Tolerance used to determine convergence of the algorithm 
#' based on the lower bound. (Default = 1e-6)
#' @param exclude Optional. An optional indicator vector of length m of 0's and 
#' 1's indicating whether to penalize a particular variable or not (0=penalize, 
#' 1=unpenalized) (Default = NULL)
#' @param add.intercept Optional. A boolean variable indicating whether or not to 
#' include an unpenalized intercept variable. (Default = TRUE)
#' @param maxit Optional. The maximum number of iterations to run the algorithm 
#' for a given solution to a penalized regression problem. (Default = 1e4)
#' @param n_orderings Optional. The number of random starts used. (Default = 10)
#' @param family Optional. The type of error model used. Currently supported modes 
#' are family="normal" and family="binomial". (Default = "normal")
#' @param scaling Optional. The type of error model used. Currently supported 
#' modes are family="normal" and family="binomial" (Default = TRUE)
#' @param return_kl Optional. A boolean variable indicating whether or not to
#' return an analysis of the null distributed features in the data-set as a 
#' function of the penalty parameter. (Default = TRUE)
#' @param estimation_type Optional. The type of estimation to perform based on 
#' the number of unique solution identified to the penalized regression problem. 
#' Valid values are estimation_type="BMA" and estimation_type="MAXIMAL" 
#' (Default = "BMA")
#' @param bma_approximation Optional. A boolean variable indicating whether to 
#' compute a full correction to the z statistic. WARNING can make the algorithm 
#' very computationally intensive for highly multi-modal posterior surfaces.
#' (Default = TRUE)
#' @param screen Optional. P-value to do marginal screening. Default is to not
#' do marginal prescreening (e.g marginal p-value of 1.0) (Default = 1.0)
#' @param post Optional. Choice of penalty parameter such that a feature will
#' have a posterior probability of 0.95 if it passes a Bonferroni correction 
#' in the multivariate model. Default is post=.95. More conservative approach 
#' would be post=0.5(Default =0.95)
#' @param already_screened Optional. If features are already screened, the 
#' marginal p-value used for screening. (Default = 1.0)
#' @param kl Optional. If features are already screened, the marginal p-value 
#' used for screening. (Default = 0.99)
#' @param l0_path Optional. The path of penalty parameters to solve the spike 
#' regression problem. If post is specified, this is computed automatically.
#' (Default =NULL)
#' @param cleanSolution Optional. This parameter determines whether a given 
#' solution is further filtered using an unpenalized model. If cleanSolution=TRUE, 
#' then the features that are significant after a Bonferroni correction given the 
#' p-values from the vbsr regression model are then tested in an unpenalized 
#' linear regression model. The p-values and z-statistics are updated using the 
#' Wald test from the unpenalized linear regression model for the features that 
#' were selected.(Default =FALSE)
#'  
#' @return A coexpression value
#' @export
sparrowZ <- function(y,x,ordering_mat=NULL, eps=1e-6, exclude=NULL, add.intercept=TRUE,
                     maxit = 1e4, n_orderings = 10, family = "normal", scaling = TRUE,
                     return_kl = TRUE, estimation_type = "BMA", bma_approximation = TRUE,
                     screen = 1.0, post=0.95, already_screened = 1.0, kl = 0.99,
                     l0_path=NULL, cleanSolution=FALSE){
  #require(vbsr)
  result <- vbsr::vbsr(y=y,X=x,ordering_mat=ordering_mat, eps=eps, exclude=exclude, 
                       add.intercept=add.intercept, maxit = maxit, 
                       n_orderings = n_orderings, family = family, scaling = scaling,
                       return_kl = return_kl, estimation_type = estimation_type, 
                       bma_approximation = bma_approximation, screen = screen,
                       post=post, already_screened = already_screened, kl = kl,
                       l0_path=l0_path, cleanSolution=cleanSolution)$z
  return(result)
}
