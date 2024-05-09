ranktransform <- function(data, thresh) rank(data)[data <= thresh]/(length(data) + 1) 
gpdtransform <- function(data, thresh, par, q) 1 - (1 - q)*pgpd(data, loc = thresh, scale = par[1], shape = par[2], lower.tail = F)

empirical_cdf <- function(data, qmarg = 0.95){ 
  u <- c()
  thresh <- quantile(data, qmarg)
  par <- gpd.fit(data, threshold = thresh, show = FALSE)$mle
  u[data <= thresh] <- ranktransform(data, thresh)
  u[data > thresh] <- gpdtransform(data[data > thresh], thresh, par, qmarg)
  return(u)
} 

#' Marginal Transformation
#' 
#' @name margtransf
#' 
#' @description
#' Marginal transformation of the random vector to standard exponential margins following \insertCite{MurphyBarltropetal2023;textual}{ReturnCurves}. 
#' 
#' @docType methods
#' 
#' @param data A matrix containing the data on the original margins.
#' @param qmarg Marginal quantile used to fit the Generalised Pareto Distribution. Default is 0.95.
#' 
#' @return A matrix containing the data on standard exponential margins.
#' 
#' @rdname marginaltransformation
#' 
#' @references \insertAllCited{}
#' 
#' @aliases margtransf
#' 
#' @examples
#' library(ReturnCurves)
#' 
#' # Generating data for illustration purposes
#' set.seed(321)
#' data <- cbind(rnorm(100), runif(100))
#' \dontrun{plot(data, pch = 20))} 
#' 
#' dataexp <- margtransf(data)
#' \dontrun{plot(dataexp, pch = 20))} 
#' 
#' @export
#' 
margtransf <- function(data, qmarg = 0.95){
  dataunif <- apply(data, 2, empirical_cdf, qmarg = qmarg)
  apply(dataunif, 2, qexp)
}  




