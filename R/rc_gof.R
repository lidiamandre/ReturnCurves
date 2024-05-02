#' Diagnostic tool for the Return Curve Estimates
#' 
#' @name rc_gof
#' 
#' @description
#' computes the goodness of fir of the return curve estimates
#' 
#' @docType methods
#' 
#' @param data matrix that contains the data, in standard exponential margins
#' @param w sequence of angles between 0 and 1; default set to a vector of 101 equally spaced angles 
#' @param method method to be used in the estimation of the angular dependence function: "hill" to use the Hill estimator, "cl" for the composite likelihood estimator
#' @param rc_origin matrix containing the return curve estimates in the original margins
#' @param q quantile to be used for the Hill estimator and/or the Heffernan and Tawn conditional extremes model; default set to 0.95
#' @param k polynomial degree for the Bernstein-Bezier polynomials used in the estimation of the angular dependence function using the composite likelihood method; default set to 7
#' @param constrained indicates whether or not to incorporate knowledge of the conditional extremes parameters; default set to "no" 
#' @param blocksize size of the blocks for the block bootstrap; default to 1 for a standard bootstrap approach
#' @param nboot number of bootstrap samples; default to 250
#' @param nangles number of angles \eqn{m} in the \eqn{(0, \pi/2)} interval; default is set to 150
#' @param alpha significance level to compute the confidence intervals
#' 
#' @return return a list with estimates for the median and lower and upper bounds of the CI for the empirical probability of lying in a survival region
#' 
#' @details to do
#'
#' @rdname rc_gof
#' 
#' @references to do
#' 
#' @aliases rc_gof
#' 
#' @examples
#' library(ReturnCurves)
#' 
#' @export
#'  
rc_gof <- function(data, w = seq(0, 1, by = 0.01), method = c("hill", "cl"), rc_origin, q = 0.95, k = 7, 
                   constrained = "no", blocksize = 1, nboot = 250, nangles = 150, alpha = 0.05){ 
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




