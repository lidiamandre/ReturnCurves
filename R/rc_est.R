#' Return Curve Estimation
#' 
#' @name rc_est
#' 
#' @description
#' computes the return curve estimation... to do
#' 
#' @docType methods
#' 
#' @param data matrix that contains the data, in standard exponential margins
#' @param w sequence of angles between 0 and 1; default set to a vector of 101 equally spaced angles 
#' @param p probability for the return curve
#' @param method method to be used in the estimation of the angular dependence function: "hill" to use the Hill estimator, "cl" for the composite likelihood estimator
#' @param q quantile to be used for the min-projection variable and Hill estimator; default set to 0.95
#' @param qalphas quantile to be used for the Heffernan and Tawn conditional extremes model; default set to 0.95
#' @param k polynomial degree for the Bernstein-Bezier polynomials used in the estimation of the angular dependence function using the composite likelihood method; default set to 7
#' @param constrained indicates whether or not to incorporate knowledge of the conditional extremes parameters; default set to "no" 
#' 
#' @return return curve estimates
#' 
#' @details This function estimates the return curve given by 
#' \deqn{RC(p):={(x, y) \in \mathbb{R}^2: \text{Pr}(X>x, Y>y)=p}.} 
#' ... talk about how it connects to the estimation of the adf and the methods used, reference the \code{\link{adf_est}} function
#' 
#' @rdname returncurve
#' 
#' @references to do
#' 
#' @aliases rc_est
#' 
#' @examples
#' library(ReturnCurves)
#' 
#' @export
#' 
rc_est <- function(data, w = seq(0, 1, by = 0.01), p, method = c("hill", "cl"), q = 0.95, qalphas = 0.95, k = 7, constrained = "no"){
  if(!method %in% c("hill", "cl")){
    stop("ADF should be estimated through the Hill estimator or Composite likelihood MLE") # write a better message here!
  }
  n <- length(w)
  xp <- qexp(1 - p)
  lambda <- adf_est(data = data, w = w, method = method, qhill = q, qalphas = qalphas, k = k, constrained = constrained)
  lambda <- properties(w, lambda)
  thresh <- sapply(w, function(i) minproj_lambda(data, i, q = q)$thresh)
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



