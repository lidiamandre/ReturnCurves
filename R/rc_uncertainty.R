#' Uncertainty of the Return Curve Estimates
#' 
#' @name rc_unc
#' 
#' @description
#' computes the uncertainty of the return curve estimates
#' 
#' @docType methods
#' 
#' @param data matrix that contains the data in the original margins
#' @param w sequence of angles between 0 and 1; default set to a vector of 101 equally spaced angles 
#' @param p probability for the return curve
#' @param method method to be used in the estimation of the angular dependence function: "hill" to use the Hill estimator, "cl" for the composite likelihood estimator
#' @param q quantile to be used for the Hill estimator and/or the Heffernan and Tawn conditional extremes model; default set to 0.95
#' @param k polynomial degree for the Bernstein-Bezier polynomials used in the estimation of the angular dependence function using the composite likelihood method; default set to 7
#' @param constrained indicates whether or not to incorporate knowledge of the conditional extremes parameters; default set to "no" 
#' @param blocksize size of the blocks for the block bootstrap; default to 1 for a standard bootstrap approach
#' @param nboot number of bootstrap samples; default to 250
#' @param nangles number of angles \eqn{m} in the \eqn{(0, \pi/2)} interval; default is set to 150
#' @param alpha significance level to compute the confidence intervals
#' 
#' @return return a list with estimates for the median, mean return curves and lower and upper bounds of the CI
#' 
#' @details to do
#'
#' @rdname rc_uncertainty
#' 
#' @references to do
#' 
#' @aliases rc_unc
#' 
#' @examples
#' library(ReturnCurves)
#' 
#' @export
#' 
rc_unc <- function(data, w = seq(0, 1, by = 0.01), p, method = c("hill", "cl"), q = 0.95, k = 7, 
                   constrained = "no", blocksize = 1, nboot = 250, nangles = 150,
                   alpha = 0.05){ 
  angles <- ((nangles:1)/(nangles + 1)) * (pi/2)
  grad <- tan(angles)
  data0 <- apply(data, 2, min)
  norms <- lapply(1:nangles, function(i) vector())
  for(i in 1:nboot){
    bootdata <- block_bootstrap_function(data = data, k = blocksize)
    bootdata_exp <- margtransf(bootdata, q = q)
    rc <- rc_est(data = bootdata_exp, w = w, p = p, method = method, q = q, k = k, constrained = constrained)
    rc_orig <- curvetransf(rc, data = bootdata, q = q)
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



