ranktransform <- function(data, thresh) rank(data)[data <= thresh]/(length(data) + 1) 
gpdtransform <- function(data, thresh, par, qmarg) 1 - (1 - qmarg)*pgpd(data, loc = thresh, scale = par[1], shape = par[2], lower.tail = F)

empirical_cdf <- function(data, qmarg = 0.95) { 
  compldata <- data[complete.cases(data)]
  u <- c()
  thresh <- quantile(compldata, qmarg)
  if (qmarg == 1) {
    stop("Threshold u too high, leading to no exceedances to fit the GPD.")
  }
  par <- gpd.fit(compldata, threshold = thresh, show = FALSE)$mle
  if (par[2] <= -1) {
    warning("MLE for the shape parameter of the GPD is < -1. \n Fitted endpoint is the maximum data point.")
  }
  if (par[2] < -0.5 && par[2] > -1) {
    warning("MLE for the shape parameter of the GPD is in (-1, -0.5). \n Non-regular MLE and a very short marginal tail is estimated.")
  }
  u[!is.na(data) & data <= thresh] <- ranktransform(data = compldata, thresh = thresh)
  u[!is.na(data) & data > thresh] <- gpdtransform(data = compldata[compldata > thresh], thresh = thresh, par = par, qmarg = qmarg)
  u[is.na(data)] <- NA
  return(u)
}

.margtransf.class <- setClass("margtransf.class", representation(data = "array",
                                                                 qmarg = "numeric",
                                                                 dataexp = "array"))

margtransf.class <- function(data, qmarg, dataexp){
  .margtransf.class(data = data,
                    qmarg = qmarg,
                    dataexp = dataexp)
}

setMethod("plot", signature = list("margtransf.class"), function(x, which = c("all", "hist", "ts", "joint")){
  which <- match.arg(which)
  if(!which %in% c("all", "hist", "ts", "joint")){
    stop("Plot type not implemented.\n 'which' should be 'hist' to plot the histograms of the data;\n 'ts' to plot the time series of each variable;\n 'joint' to plot the joint distribution;\n or 'all' for all available plots.")
  }
  df <- data.frame("X" = x@data[, 1], "Y" = x@data[, 2], "Xexp" = x@dataexp[, 1], "Yexp" = x@dataexp[, 2])
  plots <- list()
  if ("all" %in% which || "hist" %in% which) {
    origX <- df %>% ggplot(aes(x = X)) + geom_histogram(col = "darkred", fill = "red", alpha = 0.3, 
                                                        na.rm = TRUE, show.legend = FALSE, binwidth = 0.3) +
      theme_minimal() + labs(x = "X", y = "Frequency") + ggtitle(TeX("Original margin of $X$"))
    expX <- df %>% ggplot(aes(x = Xexp)) + geom_histogram(col = "darkred", fill = "red", alpha = 0.3, 
                                                          na.rm = TRUE, show.legend = FALSE, binwidth = 0.3) +
      theme_minimal() + labs(x = expression(X[exp]), y = "Frequency") + ggtitle(TeX("Marginal transformation of $X$"))
    origY <- df %>% ggplot(aes(x = Y)) + geom_histogram(col = "darkblue", fill = "blue", alpha = 0.3, 
                                                        na.rm = TRUE, show.legend = FALSE, binwidth = 0.3) +
      theme_minimal() + labs(x = "Y", y = "Frequency") + ggtitle(TeX("Original margin of $Y$"))
    expY <- df %>% ggplot(aes(x = Yexp)) + geom_histogram(col = "darkblue", fill = "blue", alpha = 0.3, 
                                                          na.rm = TRUE, show.legend = FALSE, binwidth = 0.3) +
      theme_minimal() + labs(x = expression(Y[exp]), y = "Frequency") + ggtitle(TeX("Marginal transformation of $Y$"))
    plots <- c(plots, list(origX, expX, origY, expY))
  }
  if ("all" %in% which || "ts" %in% which) {
    tsorigX <- df %>% ggplot(aes(x = 1:length(X), y = X)) + geom_line(na.rm = T) +
      theme_minimal() + labs(x = "Index", y = "X") + 
      ggtitle(TeX("Time series of $X$"))
    tsexpX <- df %>% ggplot(aes(x = 1:length(Xexp), y = Xexp)) + geom_line(na.rm = T) +
      theme_minimal() + labs(x = "Index", y = expression(X[exp])) + 
      ggtitle(TeX("Time series of $X_{exp}$"))
    tsorigY <- df %>% ggplot(aes(x = 1:length(Y), y = Y)) + geom_line(na.rm = T) +
      theme_minimal() + labs(x = "Index", y = "Y") + 
      ggtitle(TeX("Time series of $Y$"))
    tsexpY <- df %>% ggplot(aes(x = 1:length(Yexp), y = Yexp)) + geom_line(na.rm = T) +
      theme_minimal() + labs(x = "Index", y = expression(Y[exp])) + 
      ggtitle(TeX("Time series of $Y_{exp}$"))
    plots <- c(plots, list(tsorigX, tsexpX, tsorigY, tsexpY))
  }
  if ("all" %in% which || "joint" %in% which) {
    origjoint <- df %>% ggplot(aes(x = X, y = Y)) + geom_point(na.rm = TRUE) +
      theme_minimal() +
      ggtitle("Original margins")
    expjoint <- df %>% ggplot(aes(x = Xexp, y = Yexp)) + geom_point(na.rm = TRUE) +
      theme_minimal() + labs(x = expression(X[exp]), y = expression(Y[exp])) +
      ggtitle("Standard exponential margins")
    plots <- c(plots, list(origjoint, expjoint))
  }
  plot_grid(plotlist = plots, ncol = 2)
})


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
#' @return An object of S4 class \code{margtransf.class}. This object returns the arguments of the function and extra slot \code{dataexp} containing a matrix with the data on standard exponential margins. 
#' 
#' The \code{plot} function takes an object of S4 class \code{margtransf.class}, and a \code{which} argument specifying the type of plot desired (see \strong{Examples}):
#' \item{\code{"hist"}}{Plots the marginal distributions of the two variables on original and standard exponential margins.}
#' \item{\code{"ts"}}{Plots the time series of the two variables on original and standard exponential margins.}
#' \item{\code{"joint"}}{Plots the joint distribution of the two variables on original and standard exponential margins.}
#' \item{\code{"all"}}{Plots all the above mentioned plots (default).}
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
#' margdata <- margtransf(data)
#' 
#' # Plots the marginal distributions of X and Y on original vs standard exponential margins
#' plot(margdata, which = "hist") 
#' 
#' # Plots the time series of X and Y on original vs standard exponential margins
#' plot(margdata, which = "ts") 
#' 
#' # Plots the joint distribution of X and Y on original vs standard exponential margins
#' plot(margdata, which = "joint") 
#' 
#' # Plots all the available plots
#' plot(margdata, which = "all") 
#' 
#' \dontrun{
#' # To see the the S4 object's slots
#' str(margdata)
#' 
#' # To access the matrix with the data on standard exponential margins
#' margdata@@dataexp
#' }
#' 
#' @export
#' 
margtransf <- function(data, qmarg = 0.95){
  data <- as.matrix(data)
  if(is.null(dim(data)) || dim(data)[2] > 2){
    warnings("Estimation of the Return Curves and/or ADF are only implemented for a bivariate setting.")
  }
  if(qmarg < 0 | qmarg > 1){
    stop("Marginal quantile needs to be in [0, 1].")
  }
  result <- margtransf.class(data = data, qmarg = qmarg, dataexp = array())
  nas <- colSums(is.na(data))
  # if(any(nas > 0)){
  #   invisible(sapply(1:length(nas[nas > 0]), function(i){
  #     warning(paste0("There are ", nas[i], " missing values in margin X", i, ".\n These were removed."))
  #   }))
  # }
  if(any(nas > 0)){
    indnas <- which(nas > 0)
    for(i in indnas){
      warning(paste0("There are ", nas[i], " missing values in margin X", i, ".\n These were removed."))
    }
  }
  # dataunif <- apply(data, 2, empirical_cdf, qmarg = qmarg)
  dataunif <- matrix(NA, ncol = 2, nrow = dim(data)[1])
  for(i in 1:2){
    dataunif[, i] <- empirical_cdf(data[, i], qmarg = qmarg)
  }
  result@dataexp <- apply(dataunif, 2, qexp)
  return(result)
}  




