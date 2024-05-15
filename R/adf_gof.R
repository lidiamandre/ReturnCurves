.adf_gof.class <- setClass("adf_gof.class", representation(data = "array",
                                                           w_ind = "numeric",
                                                           w = "numeric",
                                                           lambda = "numeric",
                                                           q = "numeric",
                                                           blocksize = "numeric",
                                                           nboot = "numeric",
                                                           alpha = "numeric",
                                                           gof = "list"))

adf_gof.class <- function(data, w_ind, w, lambda, q, blocksize, nboot, alpha, gof){
  .adf_gof.class(data = data,
                 w_ind = w_ind,
                 w = w,
                 lambda = lambda,
                 q = q,
                 blocksize = blocksize,
                 nboot = nboot,
                 alpha = alpha,
                 gof = gof)
}

#' Goodness of fit of the Angular Dependence function estimates
#' 
#' @name adf_gof
#' 
#' @description \loadmathjax{}
#' Assessment of the goodness of fit of the angular dependence function estimates \mjeqn{\lambda(\omega)}{} following the procedure of \insertCite{MurphyBarltropetal2024;textual}{ReturnCurves}.
#' 
#' @docType methods
#' 
#' @param data A matrix containing the data on standard exponential margins.
#' @param w_ind Index of the ray to be considered on the goodness of fit assessment.
#' @param w Sequence of angles between \code{0} and \code{1}. Default is \code{seq(0, 1, by = 0.01)}.
#' @param lambda \loadmathjax{} Vector containing the estimates of the angular dependence function \mjeqn{\lambda(\omega)}{}.
#' @param q \loadmathjax{} Marginal quantile to be used for the min-projection variable \mjeqn{T^1}{} at angle \mjeqn{\omega}{} (see \strong{Details}). Default is \code{0.95}.
#' @param blocksize Size of the blocks for the block bootstrap procedure. If \code{1} (default), then a standard bootstrap approach is applied.
#' @param nboot Number of bootstrap samples to be taken. Default is \code{250} samples.
#' @param alpha \loadmathjax{}Significance level to compute the \mjeqn{(1-\alpha)}{} confidence intervals. Default is \code{0.05}.
#' 
#' @return Returns an object of S4 class of type adf_gof.class. 
#' 
#' @slot gof A list containing
#' \item{model}{A vector containing the model quantiles.} 
#' \item{empirical}{A vector containing the empirical quantiles.}
#' \item{lower}{A vector containing the lower bound of the confidence interval.}
#' \item{upper}{A vector containing the upper bound of the confidence interval.}
#' 
#' @details \loadmathjax{} Define the min-projection variable as \mjeqn{t^1_\omega = t_\omega - u_\omega | t_\omega > u_\omega}{}, then
#' variable \mjeqn{\lambda(\omega)T^1_\omega \sim Exp(1)}{} as \mjeqn{u_\omega \to \infty}{} for all \mjeqn{\omega \in [0,1]}{}. 
#' A good fit is shown by agreement of model and empirical quantiles, i.e. points should lie close to the line \mjeqn{y=x}{} (lie within the \mjeqn{(1-\alpha)}{} confidence band).
#' 
#' @rdname adf_gof
#' 
#' @references \insertAllCited{}
#' 
#' @aliases adf_gof
#' 
#' @examples
#' library(ReturnCurves)
#' 
#' # Generating data for illustration purposes
#' set.seed(321)
#' data <- cbind(rnorm(1000), rnorm(1000))
#' 
#' dataexp <- margtransf(data)
#' 
#' adf <- adf_est(data = dataexp, method = "hill")
#' 
#' w_ind <- 31
#' gof <- adf_gof(data = dataexp, w_ind = w_ind, lambda = adf)
#' 
#' \dontrun{
#' plot(gof$model, gof$empirical, ylab = "Empirical", xlab = "Model")
#' polygon(c(rev(gof$model), gof$model), c(rev(gof$lower), gof$upper), col = 'grey', border = NA)
#' points(gof$model, gof$empirical, pch = 20, col = "black")
#' abline(0, 1, col = 2,lwd = 3) 
#' }
#' 
#' @export
#'  
adf_gof <- function(data, w_ind, w = seq(0, 1, by = 0.01), lambda, q = 0.95,
                    blocksize = 1, nboot = 250, alpha = 0.05){
  if(w_ind > length(w)) stop("Angle not considered") # future me - change this
  if(length(lambda) != length(w)) stop("Number of angles and values estimated for the adf differ") # future me - change this
  result <- adf_gof.class(data = data, w_ind = w_ind, w = w, lambda = lambda, q = q, blocksize = blocksize,
                          nboot = nboot, alpha = alpha, gof = list())
  min_proj <- ReturnCurves:::minproj_lambda(data = data, w = w[w_ind], q_minproj = q)
  excdata <- (min_proj$minproj - min_proj$thresh)[min_proj$minproj > min_proj$thresh]
  excdata <- lambda[w_ind] * excdata
  nexcdata <- length(excdata)
  empirical_quantile <- sort(excdata)
  model_quantile <- qexp((1:nexcdata) / (nexcdata + 1), rate = 1)
  empirical_quantile_boot <- matrix(NA, nrow = nboot, ncol = length(empirical_quantile))
  for(i in 1:nboot){
    bdata <- ReturnCurves:::block_bootstrap_function(data = excdata, k = blocksize, n = nexcdata)
    empirical_quantile_boot[i, ] <- sort(bdata)
  }
  ub <- apply(empirical_quantile_boot, 2, quantile, probs = 1 - alpha/2)
  lb <- apply(empirical_quantile_boot, 2, quantile, probs = alpha/2)
  result@gof <- list("model" = model_quantile, "empirical" = empirical_quantile, "lower" = lb, "upper" = ub)
  return(result)
}


