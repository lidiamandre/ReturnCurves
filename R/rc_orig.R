curve_inverse_transform <- function(curveunif, data, qmarg = 0.95){
  # compldata <- data[complete.cases(data)]
  thresh <- quantile(data, qmarg)
  if(qmarg == 1){
    stop("Threshold u too high, leading to no exceedances to fit the GPD.")
  }
  par <- gpd.fit(data, threshold = thresh, show = FALSE)$mle
  if(par[2] <= -1){
    warning("MLE for the shape parameter of the GPD is < -1. \n Fitted endpoint is the maximum data point.")
  }
  if(par[2] < -0.5 && par[2] > -1){
    warning("MLE for the shape parameter of the GPD is in (-1, -0.5). \n Non-regular MLE and a very short marginal tail is estimated.")
  }
  nvec <- c()
  nvec[curveunif > qmarg] <- qgpd((curveunif[curveunif > qmarg] - qmarg)/(1 - qmarg), loc = thresh, scale = par[1], shape = par[2])
  nvec[curveunif <= qmarg] <- quantile(data, curveunif[curveunif <= qmarg])
  return(nvec)
}

.rc_est.class <- setClass("rc_est.class", representation(data = "array",
                                                         qmarg = "numeric",
                                                         w = "numeric",
                                                         p = "numeric",
                                                         method = "character",
                                                         q = "numeric",
                                                         qalphas = "numeric",
                                                         k = "numeric",
                                                         constrained = "logical",
                                                         tol = "numeric",
                                                         par_init = "numeric",
                                                         rc = "array"))

rc_est.class <- function(data, qmarg, w, p, method, q, qalphas, k, constrained, tol, par_init, rc){
  .rc_est.class(data = data,
                qmarg = qmarg,
                w = w,
                p = p,
                method = method,
                q = q,
                qalphas = qalphas,
                k = k,
                constrained = constrained,
                tol = tol,
                par_init = par_init,
                rc = rc)
}

setMethod("plot", signature = list("rc_est.class"), function(x){
  df <- data.frame("X" = x@data[, 1], "Y" = x@data[, 2])
  rcdf <- data.frame("rcX" = x@rc[, 1], "rcY" = x@rc[, 2])
  df %>% ggplot(aes(x = X, y = Y)) + geom_point(na.rm = T) +
    geom_line(data = rcdf, aes(x = rcX, y = rcY), col = "red", linewidth = 1) +
    theme_minimal() +
    ggtitle(TeX("Estimation of $\\hat{RC}(p)$"))
})

#' Return Curve estimation
#' 
#' @name rc_est
#' 
#' @description
#' \loadmathjax{} Estimation of the \mjeqn{p}{p}-probability return curve following \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}.
#'  
#' @docType methods
#' 
#' @param data A matrix containing the data on the original margins.
#' @param qmarg Marginal quantile to be used for the fit of the Generalised Pareto Distribution. Default is \code{0.95}.
#' @param p \loadmathjax{} Curve survival probability. Must be \mjeqn{p < 1-q}{p < 1-q} and \mjeqn{p < 1-q_\alpha}{p < 1-qalphas}.
#' @inheritParams adf_est
#' 
#' @return An object of S4 class \code{rc_est.class}. This object returns the arguments of the function and extra slot \code{rc} containing a matrix with the estimates of the Return Curve.
#' 
#' @details \loadmathjax{} Let \mjeqn{X, Y\sim Exp(1)}{}. Given a probability \mjeqn{p}{p} and a joint survival function \mjeqn{Pr(X>x, Y>y)}{}, 
#' the \mjeqn{p}{p}-probability return curve is defined as 
#' \mjdeqn{RC(p):=\left\lbrace(x, y) \in R^2: Pr(X>x, Y>y)=p\right\rbrace.}{} 
#' 
#' \mjeqn{Pr(X>x, Y>y)}{} is estimated using the angular dependence function \mjeqn{\lambda(\omega)}{} introduced by \insertCite{WadsworthTawn2013;textual}{ReturnCurves}. More details on how to estimate \mjeqn{\lambda(\omega)}{} can be found in \code{\link{adf_est}}.
#' 
#' The return curve estimation \mjeqn{\hat{RC}(p)}{} is done in standard exponential margins and then back transformed onto the original margins.
#' 
#' @note The parameter \code{qmarg} should be the same as the one used in \code{\link{margtransf}}.
#' 
#' @rdname returncurve
#' 
#' @references \insertAllCited{}
#' 
#' @aliases rc_est
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
#' prob <- 1/n
#' 
#' rc_orig <- rc_est(data = data, p = prob, method = "hill")
#' 
#' plot(rc_orig)
#' 
#' \dontrun{
#' # To see the the S4 object's slots
#' str(rc_orig)
#' 
#' # To access the matrix with the data on standard exponential margins
#' rc_orig@@rc
#' }
#' 
#' @export
rc_est <- function(data, qmarg = 0.95, w = seq(0, 1, by = 0.01), p, method = c("hill", "cl"), q = 0.95, qalphas = 0.95, k = 7, constrained = FALSE, tol = 0.001, par_init = rep(0, k - 1)){
  if(is.null(dim(data)) || dim(data)[2] > 2){
    stop("Estimation of the Return Curve is only implemented for a bivariate setting.")
  }
  if(qmarg < 0 | qmarg > 1){
    stop("Marginal quantiles need to be in [0, 1].")
  }
  if(p < 0 | p > 1){
    stop("Probability needs to be in [0, 1].")
  }
  if(p > 1 - qmarg | p > 1 - q | p > 1 - qalphas){
    warning("The curve survival probability p should not be too extreme and within the range of the data, i.e. smaller than the marginal quantiles.")
  }
  dataexp <- margtransf(data = data, qmarg = qmarg)@dataexp
  result <- rc_est.class(data = data, qmarg = qmarg, w = w, p = p, method = method, q = q, qalphas = qalphas, k = k, constrained = constrained, tol = tol, par_init = par_init, rc = array())
  rc_data <- rc_exp(data = dataexp, w = w, p = p, method = method, q_minproj = q, qalphas = qalphas, k = k, constrained = constrained, tol = tol, par_init = par_init)
  curveunif <- apply(rc_data, 2, pexp)
  data <- data[complete.cases(data), ]
  result@data <- data
  result@rc <- sapply(1:dim(curveunif)[2], function(i) curve_inverse_transform(curveunif[, i], data = data[, i], qmarg = qmarg))
  return(result)
}




