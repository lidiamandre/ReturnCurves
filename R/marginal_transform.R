ranktransform <- function(data, thresh) rank(data)[data <= thresh]/(length(data) + 1) 
gpdtransform <- function(data, thresh, par, q) 1 - (1 - q)*pgpd(data, loc = thresh, scale = par[1], shape = par[2], lower.tail = F)

empirical_cdf <- function(data, q = 0.95){ 
  u <- c()
  thresh <- quantile(data, q)
  par <- gpd.fit(data, threshold = thresh, show = FALSE)$mle
  u[data <= thresh] <- ranktransform(data, thresh)
  u[data > thresh] <- gpdtransform(data[data > thresh], thresh, par, q)
  return(u)
} 

#' Marginal Transformation
#' 
#' @name margtransf
#' 
#' @description
#' computes the marginal transformation to standard exponential through the probability integral transform... to do
#' 
#' @docType methods
#' 
#' @param data matrix that contains the data, in the original margins
#' @param q quantile to be used for the GPD; default set to 0.95
#' 
#' @return matrix containing the variables in standard exponential margins
#' 
#' @rdname marginaltransformation
#' 
#' @aliases margtransf
#' 
#' @examples
#' library(ReturnCurves)
#' 
#' @export
#' 
margtransf <- function(data, q = 0.95){
  dataunif <- apply(data, 2, empirical_cdf, q = q)
  aplly(dataunif, 2, qexp)
}  




