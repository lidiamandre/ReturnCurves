curve_inverse_transform <- function(curveunif, data, qmarg = 0.95){
  thresh <- quantile(data, qmarg)
  par <- gpd.fit(data, threshold = thresh, show = FALSE)$mle
  nvec <- c()
  nvec[curveunif > qmarg] <- qgpd((curveunif[curveunif > qmarg] - qmarg)/(1 - qmarg), loc = thresh, scale = par[1], shape = par[2])
  nvec[curveunif <= qmarg] <- quantile(data, curveunif[curveunif <= qmarg])
  return(nvec)
}

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
#' @return A matrix containing the estimates of the Return Curve.
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
#' \dontrun{
#' plot(data, pch = 20, main = "Return Curve on the original margins")
#' lines(rc_orig, col = 2, lwd = 2)
#' }
#' 
#' @export
rc_est <- function(data, qmarg = 0.95, w = seq(0, 1, by = 0.01), p, method = c("hill", "cl"), q = 0.95, qalphas = 0.95, k = 7, constrained = FALSE, tol = 0.001){
  dataexp <- margtransf(data = data, qmarg = qmarg)
  rc_data <- rc_exp(data = dataexp, w = w, p = p, method = method, q_minproj = q, qalphas = qalphas, k = k, constrained = constrained, tol = tol)
  curveunif <- apply(rc_data, 2, pexp)
  sapply(1:dim(curveunif)[2], function(i) curve_inverse_transform(curveunif[, i], data = data[, i], qmarg = qmarg))
}
