#' Goodness of fit of the Return Curve estimates
#' 
#' @name rc_gof
#' 
#' @description
#' Assessment of the goodness-of-fit of the return curve estimates following the approach of \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}.
#' 
#' @docType methods
#' 
#' @param data A matrix containing the data on the original margins.
#' @param w Sequence of angles between 0 and 1. Default is \code{seq(0, 1, by = 0.01)}. 
#' @param rc_origin A matrix containing the estimates of the return curve, on the original margins. Should be an object of function \code{\link{curvetransform}}.
#' @param blocksize Size of the blocks for the block bootstrap procedure. If 1 (default), then a standard bootstrap approach is applied.
#' @param nboot Number of bootstrap samples to be taken. Default is 250 samples.
#' @param nangles \loadmathjax{} Number of angles in the interval \mjeqn{(0, \pi/2)}{} \insertCite{MurphyBarltropetal2023}{ReturnCurves}. Default is 150 angles.
#' @param alpha \loadmathjax{} Significance level to compute the \mjeqn{(1-\alpha)}{} confidence intervals. Default is 0.05.
#' 
#' @return Returns a list containing:
#' \item{median}{A vector containing the median of the empirical probability of lying in a survival region.} 
#' \item{lower}{A vector containing the lower bound of the confidence interval.}
#' \item{upper}{A vector containing the upper bound of the confidence interval.}
#' 
#' @details \loadmathjax{} Given a return curve RC(\mjeqn{p}{p}), the probability of lying on a survival region is \mjeqn{p}{p}. 
#' For each angle \mjeqn{\theta}{} and corresponding point in the estimated return curve \mjeqn{\lbrace \hat{x}_\theta, \hat{y}_\theta \rbrace}{}, 
#' the empirical probability \mjeqn{\hat{p}}{p} of lying in the survival region is given by the proportion of points in the region
#' \mjeqn{(\hat{x}_\theta, \infty) \times (\hat{y}_\theta, \infty)}{}. 
#' The \mjeqn{(1-\alpha)}{} confidence region is obtained via a (block) bootstrapping procedure and ideally should contain the true probability \mjeqn{p}{p}. 
#'
#' @rdname rc_gof
#' 
#' @references \insertAllCited{}
#' 
#' @aliases rc_gof
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
#' prob <- 10/(dim(data)[1])
#' 
#' rc <- rc_est(data = dataexp, p = prob, method = "hill")
#' 
#' rc_orig <- curvetransf(curvedata = rc, data = data)
#'
#' gof <- rc_gof(data = data, rc_origin = rc_orig)
#' 
#' \dontrun{
#' ang <- 1:length(gof$median)
#' plot(ang, gof$median, xlab = "Angle Index", ylab = "Probability", type = "n", ylim = c(-0.001, range(gof$upper)[2] + 0.001))
#' polygon(c(rev(ang), ang), c(rev(gof$lower), gof$upper), col = 'grey80', border = NA)
#' lines(ang, gof$median, lwd = 2)
#' lines(ang, gof$upper, lty = 'dashed', col = 'blue', lwd = 2)
#' lines(ang, gof$lower, lty = 'dashed', col = 'blue', lwd = 2)
#' lines(ang, rep(prob, length(ang)), lwd = 3, col = 2)
#' }
#' 
#' @export
#'  
rc_gof <- function(data, w = seq(0, 1, by = 0.01), rc_origin, 
                   blocksize = 1, nboot = 250, nangles = 150, alpha = 0.05){ 
  n <- dim(data)[1]
  angles <- ((nangles:1)/(nangles + 1)) * (pi/2)
  grad <- tan(angles)
  data0 <- apply(data, 2, min)
  emp_prob <- lapply(1:nangles, function(i) vector())
  curve_w <- atan((rc_origin[, 2] - data0[2])/(rc_origin[, 1] - data0[1]))
  angles_est <- matrix(NA, ncol = dim(data)[2], nrow = nangles)
  for(i in 1:nangles){
    idx <- min(which(angles[i] >= curve_w))
    data1 <- rc_origin[idx, ] - data0
    data2 <- rc_origin[idx - 1, ] - data0
    s <- (data1[1] * tan(angles[i]) - data1[2])/((data2[2] - data1[2]) - (data2[1] - data1[1]) * tan(angles[i]))
    xhat <- data1[1] + s*(data2[1] - data1[1]) + data0[1]
    yhat <- data1[2] + s*(data2[2] - data1[2]) + data0[2]
    angles_est[i, ] <- c(xhat, yhat)
  }
  for(i in 1:nboot){
    bootdata <- ReturnCurves:::block_bootstrap_function(data = data, k = blocksize, n = n)
    for(j in 1:nangles){
      emp_prob[[j]][i] <- mean(bootdata[, 1] > angles_est[j, 1] & bootdata[, 2] > angles_est[j, 2])
    }
  }
  lb <- sapply(1:nangles, function(i) quantile(emp_prob[[i]], alpha/2))
  ub <- sapply(1:nangles, function(i) quantile(emp_prob[[i]], 1 - alpha/2))
  med <- sapply(1:nangles, function(i) quantile(emp_prob[[i]], 0.5))
  return(list("median" = med, "lower" = lb, "upper" = ub))
}




