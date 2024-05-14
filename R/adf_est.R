#' Estimation of the Angular Dependence function (ADF)
#' 
#' @name adf_est
#' 
#' @description \loadmathjax{}
#' Estimation of the angular dependence function \mjeqn{\lambda(\omega)}{} introduced by \insertCite{WadsworthTawn2013;textual}{ReturnCurves}.
#' 
#' @docType methods
#' 
#' @param data A matrix containing the data on standard exponential margins.
#' @param w Sequence of angles between 0 and 1. Default is \code{seq(0, 1, by = 0.01)}.
#' @param method String that indicates which method is used for the estimation of the angular dependence function. Must either be \code{"hill"}, to use the Hill estimator \insertCite{Hill1975}{ReturnCurves}, or \code{"cl"} to use the composite maximum likelihood estimator.
#' @param q \loadmathjax{} Marginal quantile used for the min-projection variable \mjeqn{T^1}{} at angle \mjeqn{\omega}{} \mjeqn{\left(t^1_\omega = t_\omega - u_\omega | t_\omega > u_\omega\right)}{}, and/or Hill estimator \insertCite{Hill1975}{ReturnCurves}. Default is 0.95.
#' @param qalphas Marginal quantile used for the Heffernan and Tawn conditional extremes model \insertCite{HeffernanTawn2004}{ReturnCurves}. Default set to 0.95.
#' @param k Polynomial degree for the Bernstein-Bezier polynomials used for the estimation of the angular dependence function with the composite likelihood method \insertCite{MurphyBarltropetal2024}{ReturnCurves}. Default set to 7.
#' @param constrained Logical. If FALSE (default) no knowledge of the conditional extremes parameters is incorporated in the angular dependence function estimation. 
#' @param tol Convergence tolerance for the composite maximum likelihood procedure. Default set to 0.0001.
#' 
#' @return A vector containing the estimates of the angular dependence function.
#' 
#' @details \loadmathjax{} The angular dependence function \mjeqn{\lambda(\omega)}{} can be estimated through a pointwise estimator, obtained with the Hill estimator, or via a smoother approach, 
#' obtained using Composite likelihood methods. Knowledge of the conditional extremes framework introduced by \insertCite{HeffernanTawn2004;textual}{ReturnCurves} can be incorporated by setting \code{"constrained"} to \code{TRUE}.
#' For more details see \insertCite{MurphyBarltropetal2024;textual}{ReturnCurves}.
#' 
#' @note \loadmathjax{} Due to its a pointwise nature, for a better estimation of \mjeqn{\lambda(\omega)}{} 
#' it is recommended a finer grid for \mjeqn{\omega}{} (e.g. \code{w = seq(0, 1, by = 0.001)}) when \code{method = "hill"}.
#' 
#' @rdname adfestimation
#' 
#' @references \insertAllCited{}
#' 
#' @aliases adf_est
#' 
#' @examples
#' library(ReturnCurves)
#'
#' # Generating data for illustration purposes
#' set.seed(321)
#' data <- cbind(rnorm(1000), rnorm(1000))
#' 
#' dataexp <- margtransf(data)
#'
#' w <- seq(0, 1, by = 0.01)
#'
#' adf <- adf_est(data = dataexp, method = "hill")
#'
#' \dontrun{
#' plot(w, pmax(w, 1-w), type = "l", lty = 2, ylim = c(min(pmax(w, 1-w)), max(adf) + 0.1))
#' lines(w, adf, col = 2, lwd = 2)
#' }
#'
#' @export
#' 
adf_est <- function(data, w = seq(0, 1, by = 0.01), method = c("hill", "cl"), q = 0.95, qalphas = 0.95, k = 7, constrained = FALSE, tol = 0.0001){
  if(!method %in% c("hill", "cl")){
    stop("ADF should be estimated through the Hill estimator or Composite likelihood MLE") # write a better message here!
  }
  if(constrained == FALSE){
    if(method == "hill"){
      lambda_hill <- sapply(w, function(i) minproj_lambda(data, w = i, q_minproj = q)$lambdahill)
      lambda_hill <- properties(w = w, lambda = lambda_hill)
      return(lambda_hill)
    }
    else{
      a <- 0
      b <- 1
      lam_end <- c(1, 1)
      basis <- bbp(w = w, k = k, a = a, b = b)$basis
      betacl <- minfunction_mle(w = w, data = data, a = a, b = b, lam_end = lam_end, k = k, q_minproj = q, tol = tol)
      lambda_cl <- basis %*% betacl
      lambda_cl <- properties(w = w, lambda = as.vector(lambda_cl))
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
      lambda_hill[!indx] <- sapply(w[!indx], function(i) minproj_lambda(data, w = i, q_minproj = qhill)$lambdahill)
      lambda_hill <- properties(w = w, lambda = lambda_hill)
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
        basis <- bbp(w = w, k = k, a = a, b = b)$basis
        lam_end <- c(max(a, 1 - a), max(b, 1 - b))
        betacl <- minfunction_mle(w = w, data = data, a = a, b = b, lam_end = lam_end, k = k, q_minproj = q, tol = tol)
        lambda_cl[!indx] <- basis %*% betacl
        lambda_cl <- properties(w = w, lambda = as.vector(lambda_cl))
        return(lambda_cl)
      }
    }
  }
}



