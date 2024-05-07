#' Goodness of fit of the ADF estimates \eqn{\lambda(\omega)}
#' 
#' @name adf_gof
#' 
#' @description
#' Assessment of the goodness of fit of the angular dependence function estimates \eqn{\lambda(\omega)} following the procedure of \insertCite{m2024}.
#' 
#' @docType methods
#' 
#' @param data A matrix or a data frame that contains the data in standard exponential margins.
#' @param w_ind Index of the ray to be considered on the goodness of fit assessment.
#' @param w Sequence of angles between 0 and 1. Default is \code{seq(0, 1, by = 0.01)}.
#' @param lambda Vector containing the estimates of the angular dependence function \eqn{\lambda{\omega}}.
#' @param q Marginal quantile to be used for the min-projection variable \eqn{T^*} (see [Details]). Default is 0.95.
#' @param blocksize Size of the blocks for the block bootstrap procedure. If 1, then a standard bootstrap approach is applied.
#' @param nboot Number of bootstrap samples to be taken. Default is 250 samples.
#' @param alpha Significance level to compute the \eqn{(1-\alpha)} confidence intervals. Default is 0.05.
#' 
#' 
#' @return Returns a list containing: \describe{
#' \item{model}{A vector containing the model quantiles.} 
#' \item{empirical}{A vector containing the empirical quantiles.}
#' \item{lower}{A vector containing the lower bound of the confidence interval.}
#' \item{upper}{A vector containing the upper bound of the confidence interval.}
#' }
#' 
#' 
#' @details Define the min projection variable as \eqn{t*_\omega = t_\omega - u_\omega | t_\omega > u_\omega}.  
#' Variable \eqn{\lambda(\omega)T^*_\omega \sim Exp(1)} as \eqn{u_\omega \to \infty} for all \eqn{\omega \in [0,1].} 
#' Therefore, a good agreement between the model and empirical quantiles should occur and a the \eqn{y=x} between the quantiles should lie within the \eqn{(1-\alpha)} confidence band.
#' The lower and upper bounds of the confidence interval are obtained through bootstrapp of the min-projection variable.
#'
#' @rdname adf_gof
#' 
#' @references \insertRef{MurphyBarltropetal2024}{m2024}
#' 
#' @aliases adf_gof
#' 
#' @examples
#' library(ReturnCurves)
#' 
#' # Generating data for illustration purposes
#' set.seed(321)
#' data <- cbind(rnorm(100), runif(100))
#' 
#' dataexp <- margtransf(data)
#' 
#' lambda <- adf_est(data = dataexp, method = "hill")
#' 
#' w_ind <- 31
#' gof <- adf_gof(data = dataexp, w_ind = w_ind, lambda = lambda)
#' 
#' \dontrun{
#' plot(gof$model, gof$empirical, ylab = "Empirical", xlab = "Model")
#' polygon(c(rev(gof$model), gof$model), c(rev(gof$lower), gof$upper), col = 'grey', border = NA)
#' points(gof$model, gof$empirical, pch = 16, col = "black")
#' abline(0, 1, col = 2,lwd = 3) 
#' }
#' 
#' @export
#'  
adf_gof <- function(data, w_ind, w = seq(0, 1, by = 0.01), lambda, q = 0.95,
                    blocksize = 1, nboot = 250, alpha = 0.05){
  if(w_ind > length(w)) stop("Angle not considered") # future me - change this
  if(length(lambda) != length(w)) stop("Number of angles and values estimated for the adf differ") # future me - change this
  min_proj <- minproj_lambda(data = data, w = w[w_ind], q = q)
  excdata <- (min_proj$minproj - min_proj$thresh)[min_proj$minproj > min_proj$thresh]
  excdata <- lambda[w_ind] * excdata
  nexcdata <- length(excdata)
  empirical_quantile <- sort(excdata)
  model_quantile <- qexp((1:nexcdata) / (nexcdata + 1), rate = 1)
  empirical_quantile_boot <- matrix(NA, nrow = nboot, ncol = length(empirical_quantile))
  for(i in 1:nboot){
    bdata <- block_bootstrap_function(data = excdata, k = blocksize)
    empirical_quantile_boot[i, ] <- sort(bdata)
  }
  ub <- apply(empirical_quantile_boot, 2, quantile, probs = 1 - alpha/2)
  lb <- apply(empirical_quantile_boot, 2, quantile, probs = alpha/2)
  return(list("model" = model_quantile, "empirical" = empirical_quantile, "lower" = lb, "upper" = ub))
}


