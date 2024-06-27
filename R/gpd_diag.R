.marggpd.class <- setClass("marggpd.class", representation(margdata = "margtransf.class",
                                                           blocksize = "numeric",
                                                           nboot = "numeric",
                                                           alpha = "numeric",
                                                           marggpd = "list"))

#' An S4 class to represent the assessment of the Generalised Pareto Distribution Fit
#'
#' @slot margdata An S4 object of class \code{margtransf.class}. See \code{\link{margtransf}} for more details. 
#' @slot blocksize Size of the blocks for the block bootstrap procedure. If \code{1} (default), then a standard bootstrap approach is applied.
#' @slot nboot Number of bootstrap samples to be taken. Default is \code{250} samples.
#' @slot alpha Significance level to compute the \mjeqn{(1-\alpha)}{}\% confidence intervals. Default is \code{0.05}.
#' @slot gof A list containing the model and empirical exponential quantiles, and the lower and upper bound of the confidence interval.
#' 
#' @keywords internal
marggpd.class <- function(margdata, blocksize, nboot, alpha, marggpd){
  .marggpd.class(margdata = margdata,
                 blocksize = blocksize,
                 nboot = nboot,
                 alpha = alpha,
                 marggpd = marggpd)
}

#' Visualisation of the assessment of the Generalised Pareto Distribution Fit
#'
#' @description Plot method for an S4 object returned by \code{\link{marggpd}}. 
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link{marggpd}}.
#' 
#' @return A ggplot object showing the QQ-plot between the model and empirical generalised Pareto distribution quantiles.
#'
#' @rdname plotmarggpd
#'
#' @aliases plot,marggpd.class
#' 
#' @keywords internal
setMethod("plot", signature = list("marggpd.class"), function(x){
  X <- Y <- model <- empirical <- NULL # NULL them out to satisfy CRAN checks
  dfX <- data.frame("model" = x@marggpd$model[[1]], "empirical" = x@marggpd$empirical[[1]])
  ploygondfX <- data.frame("X" = c(rev(x@marggpd$model[[1]]), x@marggpd$model[[1]]),
                           "Y" = c(rev(x@marggpd$lower[[1]]), x@marggpd$upper[[1]]))
  dfY <- data.frame("model" = x@marggpd$model[[2]], "empirical" = x@marggpd$empirical[[2]])
  ploygondfY <- data.frame("X" = c(rev(x@marggpd$model[[2]]), x@marggpd$model[[2]]),
                           "Y" = c(rev(x@marggpd$lower[[2]]), x@marggpd$upper[[2]]))
  qqX <- ggplot(data = ploygondfX, aes(x = X, y = Y)) + geom_polygon(fill = "grey80", col = NA) +
    geom_point(data = dfX, mapping = aes(x = model, y = empirical)) + 
    geom_abline(col = 2, linewidth = 1) + 
    labs(x = "Model quantiles", y = "Empirical quantiles") +
    theme_minimal() +
    ggtitle("GPD fit of X")
  qqY <- ggplot(data = ploygondfY, aes(x = X, y = Y)) + geom_polygon(fill = "grey80", col = NA) +
    geom_point(data = dfY, mapping = aes(x = model, y = empirical)) + 
    geom_abline(col = 2, linewidth = 1) + 
    labs(x = "Model quantiles", y = "Empirical quantiles") +
    theme_minimal() +
    ggtitle("GPD fit of Y")
  grid.arrange(qqX, qqY)
})


#' Assessing the Generalised Pareto Distribution Fit
#' 
#' @name marggpd
#' 
#' @description \loadmathjax{}
#' Assessment of the generalised Pareto distribution fit for each margin after following the marginal transformation procedure \code{\link{margtransf}}.
#' 
#' @docType methods
#' 
#' @param margdata An S4 object of class \code{margtransf.class}. See \code{\link{margtransf}} for more details. 
#' @param blocksize Size of the blocks for the block bootstrap procedure. If \code{1} (default), then a standard bootstrap approach is applied.
#' @param nboot Number of bootstrap samples to be taken. Default is \code{250} samples.
#' @param alpha Significance level to compute the \mjeqn{(1-\alpha)}{}\% tolerance intervals. Default is \code{0.05}.
#' 
#' @return An object of S4 class \code{marggpd.class}. This object returns the arguments of the function and an extra slot \code{marggpd} which is a list containing: 
#' \item{model}{A list containing the model quantiles for each variable.} 
#' \item{empirical}{A list containing the empirical quantiles for each variable.}
#' \item{lower}{A list containing the lower bounds of the confidence interval for each variable.}
#' \item{upper}{A list containing the upper bounds of the confidence interval for each variable.}
#' 
#' @details Let \mjeqn{F^{-1}_{GPD}}{} denote the inverse of the cumulative distribution function of a variable following a Generalised Pareto Distribution (GPD) and \mjeqn{X_{(i)}}{}, \mjeqn{Y_{(i)}}{}denote the \mjeqn{i}{i}-th ordered increasing statistics, \mjeqn{i = 1, \ldots, n}{}. 
#' Function \code{plot} shows QQ plots between the model and empirical exponential quantiles for both variables, i.e. points \mjeqn{\left(F^{-1}_{GPD}\left(\frac{i}{n+1}\right), X_{(i)}\right)}{} and \mjeqn{\left(F^{-1}_{GPD}\left(\frac{i}{n+1}\right), Y_{(i)}\right)}{},
#' along with the line \mjeqn{y=x}{}. Uncertainty is obtained via a (block) bootstrap procedure and shown by the grey region on the plot.
#' A good fit is shown by agreement of model and empirical quantiles, i.e. points should lie close to the line \mjeqn{y=x}{}. 
#' In addition, line \mjeqn{y = x}{} should mainly lie within the \mjeqn{(1-\alpha)}{}\% tolerance intervals.
#' 
#' @rdname marggpd
#' 
#' @references \insertAllCited{}
#' 
#' @aliases marggpd
#' 
#' @examples
#' library(ReturnCurves)
#' 
#' data(airdata)
#' 
#' n <- dim(airdata)[1]
#' 
#' margdata <- margtransf(airdata)
#' 
#' marggpd <- marggpd(margdata = margdata)
#' 
#' plot(marggpd)
#' 
#' \dontrun{
#' # To see the the S4 object's slots
#' str(marggpd)
#' 
#' # To access the list of vectors
#' marggpd@@marggpd
#' }
#' 
#' @export
#' 
marggpd <- function(margdata, nboot = 250, blocksize = 1, alpha = 0.05){
  if(!inherits(margdata, "margtransf.class")){
    stop("The margdata argument needs to be an object of class margtransf.class.")
  }
  data <- margdata@data[complete.cases(margdata@data), ]
  thresh <- margdata@thresh
  parameters <- margdata@parameters
  result <- marggpd.class(margdata = margdata, blocksize = blocksize,
                                         nboot = nboot, alpha = alpha, marggpd = list())
  excdataX <- data[, 1][data[, 1] > thresh[1]]
  excdataY <- data[, 2][data[, 2] > thresh[2]]
  nexcX <- length(excdataX)
  nexcY <- length(excdataY)
  empquantileX <- sort(excdataX)
  empquantileY <- sort(excdataY)
  modquantileX <- qgpd((1:nexcX) / (nexcX + 1), loc = thresh[1], scale = parameters[1, 1], shape = parameters[2, 1])
  modquantileY <- qgpd((1:nexcY) / (nexcY + 1), loc = thresh[2], scale = parameters[1, 2], shape = parameters[2, 2])
  empquantilebootX <- matrix(NA, nrow = nboot, ncol = nexcX)
  empquantilebootY <- matrix(NA, nrow = nboot, ncol = nexcY)
  for(i in 1:nboot){
    bdataX <- block_bootstrap_function(data = excdataX, k = blocksize, n = nexcX)
    empquantilebootX[i, ] <- sort(bdataX)
    bdataY <- block_bootstrap_function(data = excdataY, k = blocksize, n = nexcY)
    empquantilebootY[i, ] <- sort(bdataY)
  }
  ubX <- apply(empquantilebootX, 2, quantile, probs = 1 - alpha/2)
  lbX <- apply(empquantilebootX, 2, quantile, probs = alpha/2)
  ubY <- apply(empquantilebootY, 2, quantile, probs = 1 - alpha/2)
  lbY <- apply(empquantilebootY, 2, quantile, probs = alpha/2)
  result@marggpd <- list("model" = list(modquantileX, modquantileY), 
                         "empirical" = list(empquantileX, empquantileY), 
                         "lower" = list(lbX, lbY), "upper" = list(ubX, ubY))
  return(result)
}

