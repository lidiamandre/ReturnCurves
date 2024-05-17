.rc_unc.class <- setClass("rc_unc.class", representation(data = "array",
                                                         qmarg = "numeric",
                                                         w = "numeric",
                                                         p = "numeric",
                                                         method = "character",
                                                         q = "numeric",
                                                         qalphas = "numeric",
                                                         k = "numeric",
                                                         constrained = "logical",
                                                         tol = "numeric",
                                                         blocksize = "numeric",
                                                         nboot = "numeric",
                                                         nangles = "numeric",
                                                         alpha = "numeric",
                                                         unc = "list"))

rc_unc.class <- function(data, qmarg, w, p, method, q, qalphas, k, constrained, 
                         tol, blocksize, nboot, nangles, alpha, unc){
  .rc_unc.class(data = data,
                qmarg = qmarg,
                w = w,
                p = p,
                method = method,
                q = q,
                qalphas = qalphas,
                k = k,
                constrained = constrained,
                tol = tol,
                blocksize = blocksize,
                nboot = nboot,
                nangles = nangles,
                alpha = alpha,
                unc = unc)
}

setMethod("plot", signature = list("rc_unc.class", "rc_est.class"), function(x, y, median = T, mean = T){
  object <- x
  df <- data.frame("X" = x@data[, 1], "Y" = x@data[, 2])
  rcdf <- data.frame("rcX" = y@rc[, 1], "rcY" = y@rc[, 2])
  uncdf <- data.frame("medianX" = x@unc$median[, 1], "medianY" = x@unc$median[, 2], 
                      "meanX" = x@unc$mean[, 1], "meanY" = x@unc$mean[, 2], 
                      "lowerX" = x@unc$lower[, 1], "lowerY" = x@unc$lower[, 2], 
                      "upperX" = x@unc$upper[, 1], "upperY" = x@unc$upper[, 2])
  colours <- c("Estimated RC" = "red", "Median RC" = "orange", "Mean RC" = "brown", 
               "Lower Bound" = 1, "Upper Bound" = 1)
  if(median == T && mean == T){
    df %>% ggplot(aes(x = X, y = Y)) + geom_point() +
      geom_line(data = rcdf, aes(x = rcX, y = rcY, col = names(colours)[1]), linewidth = 1) +
      geom_line(data = uncdf, aes(x = meanX, y = meanY, col = names(colours)[3]), linewidth = 1) +
      geom_line(data = uncdf, aes(x = medianX, y = medianY, col = names(colours)[2]), linewidth = 1) +
      geom_line(data = uncdf, aes(x = lowerX, y = lowerY, col = names(colours)[4]), linetype = "dashed") +
      geom_line(data = uncdf, aes(x = upperX, y = upperY, col = names(colours)[5]), linetype = "dashed") +
      scale_color_manual(values = colours, 
                         guide = guide_legend(override.aes = list(linetype = c("solid", "dashed", "solid", 
                                                                               "solid", "dashed"),
                                                                  linewidth = c(1, 0.5, 1, 1, 0.5)))) +
      theme_minimal() + theme(legend.title = element_blank()) +
      ggtitle(TeX("Uncertainty of $\\hat{RC}(p)$"))
  }
  else if(median == T && mean == F){
    df %>% ggplot(aes(x = X, y = Y)) + geom_point() +
      geom_line(data = rcdf, aes(x = rcX, y = rcY, col = names(colours)[1]), linewidth = 1) +
      geom_line(data = uncdf, aes(x = medianX, y = medianY, col = names(colours)[2]), linewidth = 1) +
      geom_line(data = uncdf, aes(x = lowerX, y = lowerY, col = names(colours)[4]), linetype = "dashed") +
      geom_line(data = uncdf, aes(x = upperX, y = upperY, col = names(colours)[5]), linetype = "dashed") +
      scale_color_manual(values = colours, 
                         guide = guide_legend(override.aes = list(linetype = c("solid", "dashed", 
                                                                               "solid", "dashed"),
                                                                  linewidth = c(1, 0.5, 1, 0.5)))) +
      theme_minimal() + theme(legend.title = element_blank()) +
      ggtitle(TeX("Uncertainty of $\\hat{RC}(p)$"))
  }
  else if(median == F && mean == T){
    df %>% ggplot(aes(x = X, y = Y)) + geom_point() +
      geom_line(data = rcdf, aes(x = rcX, y = rcY, col = names(colours)[1]), linewidth = 1) +
      geom_line(data = uncdf, aes(x = meanX, y = meanY, col = names(colours)[3]), linewidth = 1) +
      geom_line(data = uncdf, aes(x = lowerX, y = lowerY, col = names(colours)[4]), linetype = "dashed") +
      geom_line(data = uncdf, aes(x = upperX, y = upperY, col = names(colours)[5]), linetype = "dashed") +
      scale_color_manual(values = colours, 
                         guide = guide_legend(override.aes = list(linetype = c("solid", "dashed", 
                                                                               "solid", "dashed"),
                                                                  linewidth = c(1, 0.5, 1, 0.5)))) +
      theme_minimal() + theme(legend.title = element_blank()) +
      ggtitle(TeX("Uncertainty of $\\hat{RC}(p)$"))
  }
  else if(median == F && mean == F){
    df %>% ggplot(aes(x = X, y = Y)) + geom_point() +
      geom_line(data = rcdf, aes(x = rcX, y = rcY, col = names(colours)[1]), linewidth = 1) +
      geom_line(data = uncdf, aes(x = lowerX, y = lowerY, col = names(colours)[4]), linetype = "dashed") +
      geom_line(data = uncdf, aes(x = upperX, y = upperY, col = names(colours)[5]), linetype = "dashed") +
      scale_color_manual(values = colours, 
                         guide = guide_legend(override.aes = list(linetype = c("solid", "dashed", "dashed"),
                                                                  linewidth = c(1, 0.5, 0.5)))) +
      theme_minimal() + theme(legend.title = element_blank()) +
      ggtitle(TeX("Uncertainty of $\\hat{RC}(p)$"))
  } 
})

#' Uncertainty of the Return Curve estimates
#' 
#' @name rc_unc
#' 
#' @description
#' Uncertainty assessment of the return curve estimates following the procedure of \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}.
#' 
#' @docType methods
#' 
#' @inheritParams rc_est
#' @param blocksize Size of the blocks for the block bootstrap procedure. If \code{1} (default), then a standard bootstrap approach is applied.
#' @param nboot Number of bootstrap samples to be taken. Default is \code{250} samples.
#' @param nangles \loadmathjax{} Number of angles \mjeqn{m}{m} in the interval \mjeqn{(0, \pi/2)}{} \insertCite{MurphyBarltropetal2023}{ReturnCurves}. Default is \code{150} angles.
#' @param alpha \loadmathjax{} Significance level to compute the \mjeqn{(1-\alpha)}{}\% confidence intervals. Default is \code{0.05}.
#' 
#' @return An object of S4 clas \code{rc_unc.class}. This object returns the arguments of the function and an extra slot \code{unc} which is a list containing:
#' \item{median}{A vector containing the median estimates of the return curve.} 
#' \item{mean}{A vector containing the mean estimates of the return curve.} 
#' \item{lower}{A vector containing the lower bound of the confidence interval.}
#' \item{upper}{A vector containing the upper bound of the confidence interval.}
#' 
#' @details \loadmathjax{} Define a set of angles \mjdeqn{\boldsymbol{\Theta}:= \left\lbrace \frac{\pi(m+1-j)}{2(m+1)} | 1\leq j\leq m\right\rbrace}{} and \mjeqn{L_\theta:=\left\lbrace(x,y)\in R^2_+ | \tan(\theta)=y/x\right\rbrace.}{}
#' For each \mjeqn{\theta\in \boldsymbol{\Theta},}{} \mjeqn{L_\theta}{} intersects the estimated \mjeqn{RC(p)}{} exactly once, i.e., \mjeqn{\lbrace\hat{x}_\theta, \hat{y}_\theta\rbrace:= \hat{RC}(p)\cap L_\theta.}{} 
#' Uncertainty of the return curve is then quantified by the distribution of \mjeqn{\hat{d}_\theta:=\left(\hat{x}^2_\theta + \hat{y}^2_\theta\right)^{1/2}}{} via a (block) bootstrap procedure. More details can be found in \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}
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
#' n <- dim(data)[1]
#' 
#' prob <- 10/n
#' 
#' rc_orig <- rc_est(data = data, p = prob, method = "hill")
#' 
#' unc <- rc_unc(data = data, p = prob, method = "hill")
#' 
#' # Plots the estimated Return Curve, the median Return Curve and the mean Return Curve 
#' plot(unc, rc_orig) 
#' 
#' # Plots the estimated Return Curve and the median Return Curve
#' plot(unc, rc_orig, mean = F) 
#' 
#' # Plots the estimated Return Curve and the mean Return Curve 
#' plot(unc, rc_orig, median = F) 
#' 
#' # Plots the estimated Return Curve 
#' plot(unc, rc_orig, median = F, mean = F) 
#' 
#' @export
#' 
rc_unc <- function(data, w = seq(0, 1, by = 0.01), p, method = c("hill", "cl"), qmarg = 0.95, q = 0.95, 
                   qalphas = 0.95, k = 7, constrained = FALSE, blocksize = 1, 
                   nboot = 250, nangles = 150, alpha = 0.05, tol = 0.0001){ 
  if(nboot < 1 | nboot %% 1 != 0){
    stop("The number of bootstrap samples needs to be a positive integer.")
  }
  if(nangles < 1 | nangles %% 1 != 0){
    stop("The number of angles needs to be a positive integer.")
  }
  if(alpha < 0 | alpha > 1){
    stop("The significance level needs to be in [0, 1].")
  }
  if(alpha > 0.5){
    warning("This will lead to a confidence interval smaller than 50%. Perhaps you mean 1-alpha.")
  }
  result <- rc_unc.class(data = data, qmarg = qmarg, w = w, p = p, method = method, 
                         q = q, qalphas = qalphas, k = k, constrained = constrained, 
                         tol = tol, blocksize = blocksize, nboot = nboot, nangles = nangles, 
                         alpha = alpha, unc = list())
  n <- dim(data)[1]
  angles <- ((nangles:1)/(nangles + 1)) * (pi/2)
  grad <- tan(angles)
  data0 <- apply(data, 2, min)
  norms <- lapply(1:nangles, function(i) vector())
  for(i in 1:nboot){
    bootdata <- ReturnCurves:::block_bootstrap_function(data = data, k = blocksize, n = n)
    rc_orig <- rc_est(data = bootdata, qmarg = qmarg, w = w, p = p, method = method, q = q, qalphas = qalphas, k = k, constrained = constrained, tol = tol)@rc
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
  result@unc <- list("median" = rc_median, "mean" = rc_mean, "lower" = rc_lb, "upper" = rc_ub)
  return(result)
}



