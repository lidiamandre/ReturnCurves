rc_est <- function(data, w, p, method, q_minproj, qalphas, k, constrained, tol){
  if(!method %in% c("hill", "cl")){
    stop("ADF should be estimated through the Hill estimator or Composite likelihood MLE") # write a better message here!
  }
  n <- length(w)
  xp <- qexp(1 - p)
  lambda <- adf_est(data = data, w = w, method = method, q = q_minproj, qalphas = qalphas, k = k, constrained = constrained, tol = tol)
  thresh <- sapply(w, function(i) minproj_lambda(data, i, q_minproj = q_minproj)$thresh)
  r <- sapply(1:n, function(i) thresh[i] - log(p/(1 - q))/lambda[i])
  x <- sapply(1:n, function(i) r[i] * w[i])
  y <- sapply(1:n, function(i) r[i] * (1 - w[i]))
  for(i in 1:length(x)){
    if(x[i] > xp){
      x[i] <- xp
    }
    if(x[i] < 0){
      x[i] <- 0
    }
    if(y[i] > xp){
      y[i] <- xp
    }
    if(y[i] < 0){
      y[i] <- 0
    }
  }
  x[1] <- 0
  y[1] <- xp
  x[n] <- xp
  y[n] <- 0
  for(i in length(w[w < 0.5]):1){
    if(x[i] > x[i + 1]){
      x[i] <- x[i + 1]
    }
    if(y[i] < y[i + 1]){
      y[i] <- y[i + 1]
    }
  }
  for(i in length(w[w > 0.5]):(length(w))){
    if(x[i] < x[i - 1]){
      x[i] <- x[i - 1]
    }
    if(y[i] > y[i - 1]){
      y[i] <- y[i - 1]
    }
  }
  return(cbind(x, y))
}

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
#' @name rc_o
#' 
#' @description
#' \loadmathjax{} Estimation of the \mjeqn{p}{p}-probability return curve following \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}.
#'  
#' @docType methods
#' 
#' @param data A matrix containing the data on the original margins.
#' @param qmarg Marginal quantile to be used for the fit of the Generalised Pareto Distribution. Default is \code{0.95}.
#' @param w Sequence of angles between \code{0} and \code{1}. Default is \code{seq(0, 1, by = 0.01)}.
#' @param p \loadmathjax{} Curve survival probability. Must be \mjeqn{p < 1-q}{p < 1-q} and \mjeqn{p < 1-q_\alpha}{p < 1-qalphas}.
#' @param method String that indicates which method is used for the estimation of the angular dependence function. Must either be \code{"hill"}, to use the Hill estimator \insertCite{Hill1975}{ReturnCurves}, or \code{"cl"} to use the composite maximum likelihood estimator. More details can be found in \code{\link{adf_est}}.
#' @param q \loadmathjax{} Marginal quantile used for the min-projection variable \mjeqn{T^1}{} at angle \mjeqn{\omega}{} \mjeqn{\left(t^1_\omega = t_\omega - u_\omega | t_\omega > u_\omega\right)}{}, and/or Hill estimator \insertCite{Hill1975}{ReturnCurves}. Default is \code{0.95}.
#' @param qalphas Marginal quantile used for the Heffernan and Tawn conditional extremes model \insertCite{HeffernanTawn2004}{ReturnCurves}. Default set to \code{0.95}.
#' @param k Polynomial degree for the Bernstein-Bezier polynomials used for the estimation of the angular dependence function with the composite likelihood method \insertCite{MurphyBarltropetal2023}{ReturnCurves}. Default set to \code{7}.
#' @param constrained Logical. If \code{FALSE} (default) no knowledge of the conditional extremes parameters is incorporated in the angular dependence function estimation. 
#' @param tol Convergence tolerance for the composite maximum likelihood procedure. Default set to \code{0.0001}.
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
#' @rdname returncurve_o
#' 
#' @references \insertAllCited{}
#' 
#' @aliases rc_o
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
#' rc_orig <- rc_o(data = data, p = prob, method = "hill")
#' 
#' \dontrun{
#' plot(data, pch = 20, main = "Return Curve on the original margins")
#' lines(rc_orig, col = 2, lwd = 2)
#' }
#' 
#' @export
rc_o <- function(data, qmarg = 0.95, w = seq(0, 1, by = 0.01), p, method = c("hill", "cl"), q = 0.95, qalphas = 0.95, k = 7, constrained = FALSE, tol = 0.001){
  dataexp <- margtransf(data = data, qmarg = qmarg)
  rc_data <- rc_est(data = dataexp, w = w, p = p, method = method, q_minproj = q, qalphas = qalphas, k = k, constrained = constrained, tol = tol)
  curveunif <- apply(rc_data, 2, pexp)
  sapply(1:dim(curveunif)[2], function(i) curve_inverse_transform(curveunif[, i], data = data[, i], qmarg = qmarg))
}
