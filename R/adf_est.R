#' Estimation of the Angular Dependence function \eqn{\lambda(\omega)}
#' 
#' @name adf_est
#' 
#' @description
#' estimates the angular dependence function \eqn{\lambda(\omega)} ... to do
#' 
#' @docType methods
#' 
#' @param data matrix that contains the data, in standard exponential margins
#' @param w sequence of angles between 0 and 1; default set to a vector of 101 equally spaced angles 
#' @param method method to be used in the estimation of the angular dependence function: "hill" to use the Hill estimator, "cl" for the composite likelihood estimator
#' @param qhill quantile to be used for the Hill estimator; default set to 0.95
#' @param qalphas quantile to be used for the Heffernan and Tawn conditional extremes model; default set to 0.95
#' @param k polynomial degree for the Bernstein-Bezier polynomials used in the estimation of the angular dependence function using the composite likelihood method; default set to 7
#' @param constrained indicates whether or not to incorporate knowledge of the conditional extremes parameters; default set to "no" 
#' 
#' @return adf estimates
#' 
#' @details This function returns the estimation of the Angular Dependence Function (ADF) denoted by \eqn{\lambda(\omega)} given in ...
#' 
#' This can be estimated through a pointwise estimator, obtained with the Hill estimator, or via a smooth approach, 
#' obtained using Composite likelihood methods. Knowledge of the conditional extremes framework introduced by Heffernan 
#' and Tawn can be incorporated by setting the argument "contrained" to "no"
#' 
#' @rdname adfestimation
#' 
#' @references to do
#' 
#' @aliases adf_est
#' 
#' @examples
#' library(ReturnCurves)
#' 
#' @export
#' 
adf_est <- function(data, w = seq(0, 1, by = 0.01), method = c("hill", "cl"), qhill = 0.95, qalphas = 0.95, k = 7, constrained = "no"){
  if(!method %in% c("hill", "cl")){
    stop("ADF should be estimated through the Hill estimator or Composite likelihood MLE") # write a better message here!
  }
  if(constrained == "no"){
    if(method == "hill"){
      lambda_hill <- sapply(w, function(i) minproj_lambda(data, w = i, q = qhill)$lambdahill)
      lambda_hill <- properties(w, lambda_hill)
      return(lambda_hill)
    }
    else{
      basis <- bbp(w = w)$basis
      betacl <- minfunction_mle(w = w, data = data)
      lambda_cl <- basis %*% betacl
      lambda_cl <- properties(w, as.vector(lambda_cl))
      return(lambda_cl)
    }
  }
  else{
    alphas <- heff_tawn_alphas(data = data, q = qalphas)
    a <- alphas[1]/(1 + alphas[1])
    b <- 1/(1 + alphas[2])
    indx <- w < a | w > b
    if(method == "hill"){
      lambda_hill <- c()
      if(sum(!indx) < 2){
        lambda_hill <- pmax(w, 1 - w)
        return(lambda_hill)
      }
      lambda_hill[indx] <- pmax(w, 1 - w)[indx]
      lambda_hill[!indx] <- sapply(w[!indx], function(i) minproj_lambda(data, w = i, q = qhill)$lambdahill)
      lambda_hill <- properties(w, lambda_hill)
      return(lambda_hill)
    }
    if(method == "cl"){
      lambda_cl <- c()
      if(sum(!indx) < 2){
        lambda_cl <- pmax(w, 1 - w)
        return(lambda_cl)
      }
      else{
        lambda_cl[indx] <- pmax(w, 1 - w)[indx]
        basis <- bbp(w = w, a = a, b = b)$basis
        lam_end <- c(max(a, 1 - a), max(b, 1 - b))
        betacl <- minfunction_mle(w = w, data = data, a = a, b = b, lam_end = lam_end)
        lambda_cl[!indx] <- basis %*% betacl
        lambda_cl <- properties(w, as.vector(lambda_cl))
        return(lambda_cl)
      }
    }
  }
}



