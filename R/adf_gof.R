.adf_gof.class <- setClass("adf_gof.class", representation(adf = "adf_est.class",
                                                           w_ind = "numeric",
                                                           blocksize = "numeric",
                                                           nboot = "numeric",
                                                           alpha = "numeric",
                                                           gof = "list"))

adf_gof.class <- function(adf, w_ind, blocksize, nboot, alpha, gof){
  .adf_gof.class(adf = adf,
                 w_ind = w_ind,
                 blocksize = blocksize,
                 nboot = nboot,
                 alpha = alpha,
                 gof = gof)
}

setMethod("plot", signature = list("adf_gof.class"), function(x){
  object <- x
  df <- as.data.frame(x@gof)
  ploygondf <- data.frame("X" = c(rev(x@gof$model), x@gof$model),
                          "Y" = c(rev(x@gof$lower), x@gof$upper))
  ploygondf %>% ggplot(aes(x = X, y = Y)) + geom_polygon(fill = "grey80", col = NA) +
    geom_point(data = df, mapping = aes(x = model, y = empirical)) + 
    geom_abline(col = 2, linewidth = 1) + 
    labs(x = "Model quantiles", y = "Empirical quantiles") +
    theme_minimal() +
    ggtitle(TeX("Goodness of fit of $\\hat{\\lambda}(\\omega)$"))
})

#' Goodness of fit of the Angular Dependence function estimates
#' 
#' @name adf_gof
#' 
#' @description \loadmathjax{}
#' Assessment of the goodness of fit of the angular dependence function estimates \mjeqn{\lambda(\omega)}{} following the procedure of \insertCite{MurphyBarltropetal2024;textual}{ReturnCurves}.
#' 
#' @docType methods
#' 
#' @param adf An S4 object of class \code{adf_est.class}. See \code{\link{adf_est}} for more details.
#' @param w_ind Index of the ray to be considered on the goodness of fit assessment.
#' @param blocksize Size of the blocks for the block bootstrap procedure. If \code{1} (default), then a standard bootstrap approach is applied.
#' @param nboot Number of bootstrap samples to be taken. Default is \code{250} samples.
#' @param alpha \loadmathjax{}Significance level to compute the \mjeqn{(1-\alpha)}{}\% confidence intervals. Default is \code{0.05}.
#' 
#' @return An object of S4 class \code{adf_gof.class}. This object returns the arguments of the function and an extra slot \code{gof} which is a list containing: 
#' \item{model}{A vector containing the model quantiles.} 
#' \item{empirical}{A vector containing the empirical quantiles.}
#' \item{lower}{A vector containing the lower bound of the confidence interval.}
#' \item{upper}{A vector containing the upper bound of the confidence interval.}
#' 
#' @details \loadmathjax{} Define the min-projection variable as \mjeqn{t^1_\omega = t_\omega - u_\omega | t_\omega > u_\omega}{}, then
#' variable \mjeqn{\lambda(\omega)T^1_\omega \sim Exp(1)}{} as \mjeqn{u_\omega \to \infty}{} for all \mjeqn{\omega \in [0,1]}{}. 
#' A good fit is shown by agreement of model and empirical quantiles, i.e. points should lie close to the line \mjeqn{y=x}{} (lie within the \mjeqn{(1-\alpha)}{} confidence band).
#' 
#' @rdname adf_gof
#' 
#' @references \insertAllCited{}
#' 
#' @aliases adf_gof
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
#' lambda <- adf_est(data = dataexp, method = "hill")
#' 
#' w_ind <- 31
#' gof <- adf_gof(adf = lambda, w_ind = w_ind)
#' 
#' plot(gof)
#' 
#' @export
#'  
adf_gof <- function(adf, w_ind, blocksize = 1, nboot = 250, alpha = 0.05){
  w <- adf@w
  data <- adf@data
  lambda <- adf@adf
  q <- adf@q
  if(w_ind > length(w)) stop("Angle not considered") # future me - change this
  if(length(lambda) != length(w)) stop("Number of angles and values estimated for the adf differ") # future me - change this
  result <- adf_gof.class(adf = adf, w_ind = w_ind, blocksize = blocksize,
                          nboot = nboot, alpha = alpha, gof = list())
  min_proj <- ReturnCurves:::minproj_lambda(data = data, w = w[w_ind], q_minproj = q)
  excdata <- (min_proj$minproj - min_proj$thresh)[min_proj$minproj > min_proj$thresh]
  excdata <- lambda[w_ind] * excdata
  nexcdata <- length(excdata)
  empirical_quantile <- sort(excdata)
  model_quantile <- qexp((1:nexcdata) / (nexcdata + 1), rate = 1)
  empirical_quantile_boot <- matrix(NA, nrow = nboot, ncol = length(empirical_quantile))
  for(i in 1:nboot){
    bdata <- ReturnCurves:::block_bootstrap_function(data = excdata, k = blocksize, n = nexcdata)
    empirical_quantile_boot[i, ] <- sort(bdata)
  }
  ub <- apply(empirical_quantile_boot, 2, quantile, probs = 1 - alpha/2)
  lb <- apply(empirical_quantile_boot, 2, quantile, probs = alpha/2)
  result@gof <- list("model" = model_quantile, "empirical" = empirical_quantile, "lower" = lb, "upper" = ub)
  return(result)
}

#' Visualisation of goodness of fit of ADF estimates
#'
#' @name plot
#'
#' @description 
#' Plot method for an S4 object returned by \code{\link{adf_est}}.
#' 
#' @docType methods
#'
#' @param x An object of an adf_est S4 class produced by \code{\link{adf_est}}.
#' 
#' @return A ggplot object.
#' 
#' @details \loadmathjax{} The plot shows a comparison between the estimates \mjeqn{\hat{\lambda}(\omega)}{} of the ADF and its theoretical lower bound.
#' 
#' @rdname plot-methods
#'
#' @aliases plot.adf
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
#' lambda <- adf_est(data = dataexp, method = "hill")
#' 
#' plot(lambda)
#' 
#' @export
