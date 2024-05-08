curve_inverse_transform <- function(curveunif, data, q = 0.95){
  thresh <- quantile(data, q)
  par <- gpd.fit(data, threshold = thresh, show = FALSE)$mle
  nvec <- c()
  nvec[curveunif > q] <- qgpd((curveunif[curveunif > q] - q)/(1 - q), loc = thresh, scale = par[1], shape = par[2])
  nvec[curveunif <= q] <- quantile(data, curveunif[curveunif <= q])
  return(nvec)
}

#' Return Curve estimates in original margins
#' 
#' @name curvetransf
#' 
#' @description
#' Back transformation of the return curve estimates onto the original margins, following \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}.
#' 
#' @docType methods
#' 
#' @param curvedata A matrix or data frame containing an object of function \code{\link{rc_est}}.
#' @param data A matrix or data frame containing the data in the original margins.
#' @param q Marginal quantile to be used for the fit of the Generalised Pareto Distribution. Default is 0.95.
#' 
#' @return A matrix or data frame containing the estimates of the Return Curve in the original margins.
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
#' data <- cbind(rnorm(100), runif(100))
#' 
#' dataexp <- margtransf(data)
#' 
#' prob <- 0.001
#' 
#' rc <- rc_est(data = dataexp, p = prob, method = "hill")
#' 
#' rc_orig <- curvetransf(curvedata = rc, data = data)
#' 
#' \dontrun{
#' plot(data, pch = 20, main = "Return Curve in original margins")
#' lines(rc_origin, col = 2, lwd = 2)
#' }
#' 
#' @export
curvetransf <- function(curvedata, data, q = 0.95){
  curveunif <- apply(curvedata, 2, pexp)
  sapply(1:dim(curveunif)[2], function(i) curve_inverse_transform(curveunif[, i], data = data[, i], q = q))
}


