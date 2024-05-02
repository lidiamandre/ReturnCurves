#' Diagnostic tool for the ADF estimates \eqn{\hat\lambda(\omega)}
#' 
#' @name adf_gof
#' 
#' @description
#' computes the goodness of fit of the angular dependence function estimates \eqn{\hat\lambda(\omega)}
#' 
#' @docType methods
#' 
#' @param data matrix that contains the data, in standard exponential margins
#' @param w_ind index of the ray to be considered while assessing the estimates
#' @param w sequence of angles between 0 and 1; default set to a vector of 101 equally spaced angles 
#' @param lambda estimates of the angular dependence function
#' @param q quantile to be used for the Hill estimator and/or the Heffernan and Tawn conditional extremes model; default set to 0.95
#' @param blocksize size of the blocks for the block bootstrap; default to 1 for a standard bootstrap approach
#' @param nboot number of bootstrap samples; default to 250
#' @param alpha significance level to compute the confidence intervals
#' 
#' @return return a list containing the model and empirical quantiles, and the lower and upper bounds for the confidence interval
#' 
#' @details to do
#'
#' @rdname adf_gof
#' 
#' @references to do
#' 
#' @aliases adf_gof
#' 
#' @examples
#' library(ReturnCurves)
#' 
#' @export
#'  
adf_gof <- function(data, w_ind, w = seq(0, 1, by = 0.01), lambda, q = 0.95,
                    blocksize = 1, nboot = 250, alpha = 0.05){
  if(w_ind > length(w)) stop("Angle not considered") # future me - change this
  if(length(lambda) != length(w)) stop("Number of angles and values estimated for the adf differ") # future me - change this
  min_proj <- minproj_lambda(data = data, w = w[w_ind], q = q)
  excdata <- (min_proj$minproj - min_proj$thresh)[min_proj$minproj > min_proj$thresh]
  nexcdata <- length(excdata)
  empirical_quantile <- sort(excdata)
  model_quantile <- qexp((1:nexcdata) / (nexcdata + 1), rate = lambda[w_ind])
  empirical_quantile_boot <- matrix(NA, nrow = nboot, ncol = length(empirical_quantile))
  for(i in 1:nboot){
    bdata <- block_bootstrap_function(data = excdata, k = blocksize)
    empirical_quantile_boot[i, ] <- sort(bdata)
  }
  ub <- apply(empirical_quantile_boot, 2, quantile, probs = 1 - alpha/2)
  lb <- apply(empirical_quantile_boot, 2, quantile, probs = alpha/2)
  return(list("model" = model_quantile, "empirical" = empirical_quantile, "lower" = lb, "upper" = ub))
}


