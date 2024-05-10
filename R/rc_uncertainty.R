#' Uncertainty of the Return Curve estimates
#' 
#' @name rc_unc
#' 
#' @description
#' Uncertainty assessment of the return curve estimates following the procedure of \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}.
#' 
#' @docType methods
#' 
#' @param data A matrix which contains the data on the original margins.
#' @param w Sequence of angles between 0 and 1. Default is \code{seq(0, 1, by = 0.01)}. 
#' @param p Curve survival probability.
#' @param method String that indicates which method is used for the estimation of the angular dependence function. Must either be \code{"hill"}, to use the Hill estimator \insertCite{Hill1975}{ReturnCurves}, or \code{"cl"} to use the composite likelihood estimator approaches. More details can be found in \code{\link{adf_est}}.
#' @param qmarg Marginal quantile used to fit the Generalised Pareto Distribution. Default is 0.95.
#' @param q \loadmathjax{} Marginal quantile used for the min-projection variable \mjeqn{T^1}{} at angle \mjeqn{\omega}{} \mjeqn{\left(t^1_\omega = t_\omega - u_\omega | t_\omega > u_\omega\right)}{}, and/or Hill estimator \insertCite{Hill1975}{ReturnCurves}. Default is 0.95.
#' @param qalphas Marginal quantile used for the Heffernan and Tawn conditional extremes model \insertCite{HeffernanTawn2004}{ReturnCurves}. Default set to 0.95.
#' @param k Polynomial degree for the Bernstein-Bezier polynomials used for the estimation of the angular dependence function using the composite likelihood method \insertCite{MurphyBarltropetal2023}{ReturnCurves}. Default set to 7.
#' @param constrained Logical. If FALSE (default) no knowledge of the conditional extremes parameters is incorporated in the angular dependence function estimation. 
#' @param blocksize Size of the blocks for the block bootstrap procedure. If 1 (default), then a standard bootstrap approach is applied.
#' @param nboot Number of bootstrap samples to be taken. Default is 250 samples.
#' @param nangles \loadmathjax{} Number of angles \mjeqn{m}{m} in the interval \mjeqn{(0, \pi/2)}{} \insertCite{MurphyBarltropetal2023}{ReturnCurves}. Default is 150 angles.
#' @param alpha \loadmathjax{} Significance level to compute the \mjeqn{(1-\alpha)}{} confidence intervals. Default is 0.05.
#' 
#' @return Returns a list containing: 
#' \item{median}{A vector containing the median estimates of the return curve.} 
#' \item{mean}{A vector containing the mean estimates of the return curve.} 
#' \item{lower}{A vector containing the lower bound of the confidence interval.}
#' \item{upper}{A vector containing the upper bound of the confidence interval.}
#'
#' @details \loadmathjax{} Define a set of angles \mjdeqn{\boldsymbol{\Theta}:= \left\lbrace \frac{\pi(m+1-j)}{2(m+1)} | 1\leq j\leq m\right\rbrace}{} and \mjeqn{L_\theta:=\left\lbrace(x,y)\in R^2_+ | \tan(\theta)=y/x\right\rbrace.}{}
#' For each \mjeqn{\theta\in \boldsymbol{\Theta},}{} \mjeqn{L_\theta}{} intersects the estimated \mjeqn{RC(p)}{} exactly once, i.e., \mjeqn{\lbrace\hat{x}_\theta, \hat{y}_\theta\rbrace:= \hat{RC}(p)\cap L_\theta.}{} 
#' Uncertainty of the return curve is then quantified by the distribution of \mjeqn{\hat{d}_\theta:=\left(\hat{x}^2_\theta + \hat{y}^2_\theta\right)^{1/2}}{} via a (block) bootstrap procedure. More details can be found in \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}
#' 
#' 
#' 
#' @rdname rc_uncertainty
#' 
#' @references \insertAllCited{}
#' 
#' @aliases rc_unc
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
#' unc <- rc_unc(data = data, p = prob, method = "hill")
#' 
#' \dontrun{
#' plot(data, xlab = "X", ylab = "Y", pch = 20, col = "grey")
#' lines(rc_orig, lwd = 2, col = 2)
#' lines(unc$median, lwd = 2, col = "orange") # to plot median estimates
#' lines(unc$mean, lwd = 2, col = "orange") # to plot mean estimates
#' lines(unc$upper, lty = 'dashed', lwd = 2)
#' lines(unc$lower, lty = 'dashed', lwd = 2)
#' }
#' 
#' @export
#' 
rc_unc <- function(data, w = seq(0, 1, by = 0.01), p, method = c("hill", "cl"), qmarg = 0.95, q = 0.95, 
                   qalphas = 0.95, k = 7, constrained = FALSE, blocksize = 1, 
                   nboot = 250, nangles = 150, alpha = 0.05){ 
  angles <- ((nangles:1)/(nangles + 1)) * (pi/2)
  grad <- tan(angles)
  data0 <- apply(data, 2, min)
  norms <- lapply(1:nangles, function(i) vector())
  for(i in 1:nboot){
    bootdata <- ReturnCurves:::block_bootstrap_function(data = data, k = blocksize)
    bootdata_exp <- margtransf(bootdata, q = qmarg)
    rc <- rc_est(data = bootdata_exp, w = w, p = p, method = method, q = q, qalphas = qalphas, k = k, constrained = constrained)
    rc_orig <- curvetransf(rc, data = bootdata, qmarg = qmarg)
    rc_orig <- rbind(c(data0[1], rc_orig[1, 2]), rc_orig, c(rc_orig[dim(rc_orig)[1], 1], data0[2]))
    curve_w <- atan((rc_orig[, 2] - data0[2])/(rc_orig[, 1] - data0[1]))
    for(j in 1:nangles){
      idx <- min(which(angles[j] >= curve_w))
      data1 <- rc_orig[idx, ] - data0
      data2 <- rc_orig[idx - 1, ] - data0
      s <- (data1[1] * tan(angles[j]) - data1[2])/((data2[2] - data1[2]) - (data2[1] - data1[1]) * tan(angles[j]))
      xhat <- data1[1] + s*(data2[1] - data1[1])
      yhat <- data1[2] + s*(data2[2] - data1[2])
      norms[[j]][i] <- sqrt(xhat^2 + yhat^2)
    }
  }
  lb <- sapply(1:nangles, function(i) quantile(norms[[i]], alpha/2))
  ub <- sapply(1:nangles, function(i) quantile(norms[[i]], 1 - alpha/2))
  med <- sapply(1:nangles, function(i) quantile(norms[[i]], 0.5))
  mea <- sapply(1:nangles, function(i) mean(norms[[i]]))
  rc_mean <- cbind(mea/sqrt(1 + grad^2) + data0[1], grad * (mea/sqrt(1 + grad^2)) + data0[2])
  rc_median <- cbind(med/sqrt(1 + grad^2) + data0[1], grad * (med/sqrt(1 + grad^2)) + data0[2])
  rc_lb <- cbind(lb/sqrt(1 + grad^2) + data0[1], grad * (lb/sqrt(1 + grad^2)) + data0[2])
  rc_ub <- cbind(ub/sqrt(1 + grad^2) + data0[1], grad * (ub/sqrt(1 + grad^2)) + data0[2])
  return(list("median" = rc_median, "mean" = rc_mean, "lower" = rc_lb, "upper" = rc_ub))
}



