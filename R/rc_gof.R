#' Goodness of fit of the Return Curve estimates
#' 
#' @name rc_gof
#' 
#' @description
#' Assessment of the goodness-of-fit of the return curve estimates following .
#' 
#' @docType methods
#' 
#' @param data A matrix or data frame containing the data in standard exponential margins.
#' @param w Sequence of angles between 0 and 1. Default is \code{seq(0, 1, by = 0.01)}. 
#' @param rc_origin A matrix or data frame containing the estimates of the return curve, in the original margins.
#' @param blocksize Size of the blocks for the block bootstrap procedure. If 1, then a standard bootstrap approach is applied.
#' @param nboot Number of bootstrap samples to be taken. Default is 250 samples.
#' @param nangles Number of angles \eqn{m} in the \eqn{(0, \pi/2)} interval (see ). Default is 150 angles.
#' @param alpha Significance level to compute the \eqn{(1-\alpha)} confidence intervals. Default is 0.05.
#' 
#' @return Returns a list containing: \describe{
#' \item{median}{A vector containing the median of the empirical probability of lying in a survival region.} 
#' \item{lower}{A vector containing the lower bound of the confidence interval.}
#' \item{upper}{A vector containing the upper bound of the confidence interval.}
#' }
#' 
#' @details Given a return curve RC(p), the probability of lying on a survival region is p. 
#' For each angle \eqn{\theta} and corresponding point in the estimated return curve \eqn{{x_\theta, y_\theta}}, 
#' the empirical probability \eqn{\hat{p}} of lying in the survival region is given by the proportion of points in the region
#' \eqn{(x_\theta, \infty) x (y_\theta, \infty)}.
#' The true value \eqn{p} should be contained within the \eqn{(1-\alpha)} confidence region. 
#'
#' @rdname rc_gof
#' 
#' @references \insertRef{MurphyBarltropetal2023}{ReturnCurves}
#' 
#' @aliases rc_gof
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
#' prob <- 1e^-3
#' 
#' rc <- rc_est(data = dataexp, p = prob, method = "hill")
#' 
#' rc_orig <- curvetransf(curvedata = rc, data = data)
#'
#' gof <- rc_gof(data = dataexp, rc_origin = rc_orig)
#' 
#' \dontrun{
#' ang <- 1:length(gof$median)
#' plot(ang, gof$median, xlab = "Angle Index", ylab = "Probability")
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
    bootdata <- block_bootstrap_function(data = data, k = blocksize)
    for(j in 1:nangles){
      emp_prob[[j]][i] <- mean(bootdata[, 1] > angles_est[j, 1] & bootdata[, 2] > angles_est[j, 2])
    }
  }
  lb <- sapply(1:nangles, function(i) quantile(emp_prob[[i]], alpha/2))
  ub <- sapply(1:nangles, function(i) quantile(emp_prob[[i]], 1 - alpha/2))
  med <- sapply(1:nangles, function(i) quantile(emp_prob[[i]], 0.5))
  return(list("median" = med, "lower" = lb, "upper" = ub))
}




