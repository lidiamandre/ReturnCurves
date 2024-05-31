.rc_gof.class <- setClass("rc_gof.class", representation(retcurve = "rc_est.class",
                                                         blocksize = "numeric",
                                                         nboot = "numeric",
                                                         alpha = "numeric",
                                                         gof = "list"))

rc_gof.class <- function(retcurve, blocksize, nboot, alpha, gof){
  .rc_gof.class(retcurve = retcurve,
                blocksize = blocksize,
                nboot = nboot,
                alpha = alpha,
                gof = gof)
}

setMethod("plot", signature = list("rc_gof.class"), function(x){
  df <- data.frame("angles" = 1:length(x@gof$median), x@gof, "pX" = c(rev(1:length(x@gof$median)), 1:length(x@gof$median)),
                   "pY" = c(rev(x@gof$lower), x@gof$upper), "prob" = rep(x@retcurve@p, length(x@gof$median)))
  coloursl <- c("Confidence interval" = 1, "Median estimate" = 1, "True probability" = 2)
  ggplot(data = df, aes(x = pX, y = pY)) + geom_polygon(fill = "grey80", col = NA) +
    geom_line(aes(x = angles, y = median, col = names(coloursl)[2])) +
    geom_line(aes(x = angles, y = upper, col = names(coloursl)[1]), linetype = "dashed") +
    geom_line(aes(x = angles, y = lower, col = names(coloursl)[1]), linetype = "dashed") + 
    geom_line(aes(x = angles, y = prob, col = names(coloursl)[3])) +
    labs(x = "Angle Index", y = "Probability") +
    scale_color_manual(values = coloursl,
                       guide = guide_legend(override.aes = list(linetype = c("dashed", "solid", "solid"))))  +
    ylim(c(-0.001, range(df$upper)[2] + 0.001)) + 
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom") +
    ggtitle(expression("Goodness of fit of" ~ hat(RC)(p)))
  
})

#' Goodness of fit of the Return Curve estimates
#' 
#' @name rc_gof
#' 
#' @description
#' Assessment of the goodness-of-fit of the return curve estimates following the approach of \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}.
#' 
#' @docType methods
#' 
#' @param retcurve An S4 object of class \code{rc_est.class}. See \code{\link{rc_est}} for more details.
#' @inheritParams rc_unc
#'  
#' @return An object of S4 class \code{rc_gof.class}. This object returns the arguments of the function and an extra slot \code{gof} which is a list containing:
#' \item{median}{A vector containing the median of the empirical probability of lying in a survival region.} 
#' \item{lower}{A vector containing the lower bound of the confidence interval.}
#' \item{upper}{A vector containing the upper bound of the confidence interval.}
#' 
#' @details \loadmathjax{} Given a return curve RC(\mjeqn{p}{p}), the probability of lying in a survival region is \mjeqn{p}{p}. 
#' For each angle \mjeqn{\theta}{} and corresponding point in the estimated return curve \mjeqn{\lbrace \hat{x}_\theta, \hat{y}_\theta \rbrace}{}, 
#' the empirical probability \mjeqn{\hat{p}}{p} of lying in the survival region is given by the proportion of points in the region
#' \mjeqn{(\hat{x}_\theta, \infty) \times (\hat{y}_\theta, \infty)}{}. Due to the use of empirical probabilities, the value of \mjeqn{p}{p} should be within the range of the data and not too extreme. 
#' 
#' The \mjeqn{(1-\alpha)}{}\% confidence region is obtained via a (block) bootstrapping procedure and ideally should contain the true probability \mjeqn{p}{p}. 
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
#' n <- dim(data)[1]
#' 
#' prob <- 10/n
#' 
#' margdata <- margtransf(data)
#' 
#' rc_orig <- rc_est(margdata = margdata, p = prob, method = "hill")
#'
#' gof <- rc_gof(retcurve = rc_orig)
#' 
#' plot(gof)
#' 
#' \dontrun{
#' # To see the the S4 object's slots
#' str(gof)
#' 
#' # To access the list of vectors
#' gof@@gof
#' }
#' 
#' @export
#'  
rc_gof <- function(retcurve, blocksize = 1, nboot = 250, nangles = 150, alpha = 0.05){ 
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
  result <- rc_gof.class(retcurve = retcurve, blocksize = blocksize, nboot = nboot, alpha = alpha, gof = list())
  rc_origin <- result@retcurve@rc
  data <- result@retcurve@data
  w <- result@retcurve@w
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
  result@gof <- list("median" = med, "lower" = lb, "upper" = ub)
  return(result)
}




