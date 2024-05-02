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
#' @param p probability for the return curve
#' @param w sequence of angles between 0 and 1; default set to a vector of 101 equally spaced angles 
#' @param method method to be used in the estimation of the angular dependence function: "hill" to use the Hill estimator, "cl" for the composite likelihood estimator
#' @param q quantile to be used for the Hill estimator and/or the Heffernan and Tawn conditional extremes model; default set to 0.95
#' @param k polynomial degree for the Bernstein-B\'ezier polynomials used in the estimation of the angular dependence function using the composite likelihood method; default set to 7
#' @param constrained indicates whether or not to incorporate knowledge of the conditional extremes parameters; default set to "no" 
#' 
#' @return return curve estimation
#' 
#' @details \loadmathjax{} This function estimates the return curve given by 
#' \mjdeqn{RC(p):=\brace(x, y) \in \mathbb{R}^2: \text{Pr}(X>x, Y>y)=p\brace.} ... talk about how it connects to the estimation of the adf and the methods used, reference the est_lamb function
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
rc_est <- function(data, w, p, method = c("hill", "cl"), q = 0.95, k = 7, constrained = "no"){
  if(!method %in% c("hill", "cl")){
    stop("ADF should be estimated through the Hill estimator or Composite likelihood MLE") # write a better message here!
  }
  n <- length(w)
  xp <- qexp(1 - p)
  lambda <- est_lamb(data = data, w = w, method = method, q = q, k = k, constrained = constrained)
  lambda <- properties(w, lambda)
  thresh <- sapply(w, function(i) minproj_lambda(data, i)$thresh)
  r <- sapply(1:n, function(i) thresh[i] - log(p/(1-q))/lambda[i])
  x <- sapply(1:n, function(i) r[i] * w[i])
  y <- sapply(1:n, function(i) r[i] * (1 - w[i]))
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
#' @rdname rc_origin
#' 
#' @references to do
#' 
#' @aliases curvetransf
#' 
#' @examples
#' library(ReturnCurves)
#' 
curvetransf <- function(curvedata, data, q = 0.95){
  curveunif <- apply(curvedata, 2, pexp)
  sapply(1:dim(curveunif)[2], function(i) curve_inverse_transform(curveunif[, i], data = data[, i], q = q))
}


