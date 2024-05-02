curve_inverse_transform <- function(curveunif, data, q = 0.95){
  thresh <- quantile(data, q)
  par <- gpd.fit(data, threshold = thresh, show = FALSE)$mle
  nvec <- c()
  nvec[curveunif > q] <- qgpd((curveunif[curveunif > q] - q)/(1 - q), loc = thresh, scale = par[1], shape = par[2])
  nvec[curveunif <= q] <- quantile(data, curveunif[curveunif <= q])
  return(nvec)
}

#' Return Curve in Original margins
#' 
#' @name curvetransf
#' 
#' @description
#' computes the return curve estimation in the original margins following the equation in the paper
#' 
#' @docType methods
#' 
#' @param curvedata matrix containing an object of function \code{\link{rc_est}}
#' @param data matrix containing the data in the original margins
#' 
#' @return return curve estimation in original margins
#' 
#' @rdname returncurve_original
#' 
#' @references to do
#' 
#' @aliases curvetransf
#' 
#' @examples
#' library(ReturnCurves)
#' 
#' @export
curvetransf <- function(curvedata, data, q = 0.95){
  curveunif <- apply(curvedata, 2, pexp)
  sapply(1:dim(curveunif)[2], function(i) curve_inverse_transform(curveunif[, i], data = data[, i], q = q))
}


