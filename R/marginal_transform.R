ranktransform <- function(data, thresh) rank(data)[data <= thresh]/(length(data) + 1) 
gpdtransform <- function(data, thresh, par, qmarg) 1 - (1 - qmarg)*pgpd(data, loc = thresh, scale = par[1], shape = par[2], lower.tail = F)

empirical_cdf <- function(data, qmarg = 0.95){ 
  u <- c()
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
  u[data <= thresh] <- ranktransform(data = data, thresh = thresh)
  u[data > thresh] <- gpdtransform(data = data[data > thresh], thresh = thresh, par = par, qmarg = qmarg)
  return(u)
} 

#' Marginal Transformation
#' 
#' @name margtransf
#' 
#' @description
#' Marginal transformation of a random vector to standard exponential margins following \insertCite{ColesTawn1991;textual}{ReturnCurves}. 
#' 
#' @docType methods
#' 
#' @param data A matrix containing the data on the original margins.
#' @param qmarg Marginal quantile used to fit the Generalised Pareto Distribution (GPD). Default is \code{0.95}.
#' 
#' @return A matrix containing the data on standard exponential margins.
#' 
#' @details \loadmathjax{} Given a threshold value \mjeqn{u}{u}, a stationary random vector 
#' is transformed by using the empirical cumulative distribution function 
#' (cdf) below \mjeqn{u}{u}, and a GPD fit above \mjeqn{u}{u}.    
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
#' data <- cbind(rnorm(1000), rnorm(1000))
#' 
#' \dontrun{plot(data, pch = 20)} 
#' 
#' dataexp <- margtransf(data)
#' 
#' \dontrun{plot(dataexp, pch = 20)} 
#' 
#' @export
#' 
margtransf <- function(data, qmarg = 0.95){
  if(qmarg < 0 | qmarg > 1){
    stop("Marginal quantile needs to be in [0, 1].")
  }
  dataunif <- apply(data, 2, empirical_cdf, qmarg = qmarg)
  apply(dataunif, 2, qexp)
}  





