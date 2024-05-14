#' Return Curve estimation
#' 
#' @name rc_est
#' 
#' @description
#' \loadmathjax{} Estimation of the \mjeqn{p}{p}-probability return curve following \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}.
#' 
#' @docType methods
#' 
#' @param data A matrix containing the data on standard exponential margins.
#' @param w Sequence of angles between \code{0} and \code{1}. Default is \code{seq(0, 1, by = 0.01)}.
#' @param p \loadmathjax{} Curve survival probability. Must be \mjeqn{p < 1-q}{p < 1-q} and \mjeqn{p < 1-q_\alpha}{p < 1-qalphas}.
#' @param method String that indicates which method is used for the estimation of the angular dependence function. Must either be \code{"hill"}, to use the Hill estimator \insertCite{Hill1975}{ReturnCurves}, or \code{"cl"} to use the composite maximum likelihood estimator. More details can be found in \code{\link{adf_est}}.
#' @param q \loadmathjax{} Marginal quantile used for the min-projection variable \mjeqn{T^1}{} at angle \mjeqn{\omega}{} \mjeqn{\left(t^1_\omega = t_\omega - u_\omega | t_\omega > u_\omega\right)}{}, and/or Hill estimator \insertCite{Hill1975}{ReturnCurves}. Default is \code{0.95}.
#' @param qalphas Marginal quantile used for the Heffernan and Tawn conditional extremes model \insertCite{HeffernanTawn2004}{ReturnCurves}. Default set to \code{0.95}.
#' @param k Polynomial degree for the Bernstein-Bezier polynomials used for the estimation of the angular dependence function with the composite likelihood method \insertCite{MurphyBarltropetal2023}{ReturnCurves}. Default set to \code{7}.
#' @param constrained Logical. If \code{FALSE} (default) no knowledge of the conditional extremes parameters is incorporated in the angular dependence function estimation. 
#' @param tol Convergence tolerance for the composite maximum likelihood procedure. Default set to \code{0.0001}.
#' 
#' @return A matrix containing the estimates of the Return Curve on standard exponential margins.
#' 
#' @details \loadmathjax{} Given a probability \mjeqn{p}{p} and a joint survival function \mjeqn{Pr(X>x, Y>y)}{}, 
#' the \mjeqn{p}{p}-probability return curve is defined as 
#' \mjdeqn{RC(p):=\left\lbrace(x, y) \in R^2: Pr(X>x, Y>y)=p\right\rbrace.}{} 
#' 
#' \mjeqn{Pr(X>x, Y>y)}{} is estimated using the angular dependence function \mjeqn{\lambda(\omega)}{} introduced by \insertCite{WadsworthTawn2013;textual}{ReturnCurves}. More details on how to estimate \mjeqn{\lambda(\omega)}{} can be found in \code{\link{adf_est}}.
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
#' dataexp <- margtransf(data)
#' 
#' prob <- 10/n
#' 
#' rc <- rc_est(data = dataexp, p = prob, method = "hill")
#' 
#' \dontrun{
#' plot(dataexp, pch = 20, main = "Return Curve on exponential margins")
#' lines(rc, col = 2, lwd = 2)
#' }
#' 
#' @export
#' 
rc_est <- function(data, w = seq(0, 1, by = 0.01), p, method = c("hill", "cl"), q = 0.95, qalphas = 0.95, k = 7, constrained = FALSE, tol = 0.001){
  if(!method %in% c("hill", "cl")){
    stop("ADF should be estimated through the Hill estimator or Composite likelihood MLE") # write a better message here!
  }
  n <- length(w)
  xp <- qexp(1 - p)
  lambda <- adf_est(data = data, w = w, method = method, q = q, qalphas = qalphas, k = k, constrained = constrained, tol = tol)
  thresh <- sapply(w, function(i) minproj_lambda(data, i, q_minproj = q)$thresh)
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



