curve_inverse_transform <- function(curveunif, data, qmarg = 0.95){
  thresh <- quantile(data, qmarg)
  par <- gpd.fit(data, threshold = thresh, show = FALSE)$mle
  nvec <- c()
  nvec[curveunif > qmarg] <- qgpd((curveunif[curveunif > qmarg] - qmarg)/(1 - qmarg), loc = thresh, scale = par[1], shape = par[2])
  nvec[curveunif <= qmarg] <- quantile(data, curveunif[curveunif <= qmarg])
  return(nvec)
}

#' Return Curve estimates on original margins
#' 
#' @name curvetransf
#' 
#' @description
#' Back transformation of the return curve estimates onto the original margins following \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}.
#' 
#' @docType methods
#' 
#' @param curvedata A matrix containing an object of function \code{\link{rc_est}}.
#' @param data A matrix containing the data on the original margins.
#' @param qmarg Marginal quantile to be used for the fit of the Generalised Pareto Distribution. Default is \code{0.95}.
#' 
#' @return A matrix containing the estimates of the Return Curve on the original margins.
#' 
#' @note The parameter \code{qmarg} should be the same as the one used in \code{\link{margtransf}}.
#' 
#' @rdname returncurve_original
#' 
#' @references \insertAllCited{}
#' 
#' @aliases curvetransf
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
#' prob <- 1/n
#' 
#' rc <- rc_est(data = dataexp, p = prob, method = "hill")
#' 
#' rc_orig <- curvetransf(curvedata = rc, data = data)
#' 
#' \dontrun{
#' plot(data, pch = 20, main = "Return Curve on the original margins")
#' lines(rc_orig, col = 2, lwd = 2)
#' }
#' 
#' @export
curvetransf <- function(curvedata, data, qmarg = 0.95){
  curveunif <- apply(curvedata, 2, pexp)
  sapply(1:dim(curveunif)[2], function(i) curve_inverse_transform(curveunif[, i], data = data[, i], qmarg = qmarg))
}


