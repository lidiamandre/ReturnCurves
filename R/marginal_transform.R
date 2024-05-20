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

setMethod("plot", signature = list("margtransf.class"), function(x, joint = F){
  df <- data.frame("X" = x@data[, 1], "Y" = x@data[, 2], "Xexp" = x@dataexp[, 1], "Yexp" = x@dataexp[, 2])
  if(joint == F){
    origX <- df %>% ggplot(aes(x = X)) + geom_histogram(col = "darkred", fill = "red", alpha = 0.3, 
                                                        na.rm = T, show.legend = F) +
      theme_minimal() + labs(x = "X", y = "Frequency") + ggtitle(TeX("Original margin of $X$"))
    expX <- df %>% ggplot(aes(x = Xexp)) + geom_histogram(col = "darkred", fill = "red", alpha = 0.3, 
                                                          na.rm = T, show.legend = F) +
      theme_minimal() + labs(x = expression(X[exp]), y = "Frequency") + ggtitle(TeX("Marginal transformation of $X$"))
    origY <- df %>% ggplot(aes(x = Y)) + geom_histogram(col = "darkblue", fill = "blue", alpha = 0.3, 
                                                        na.rm = T, show.legend = F) +
      theme_minimal() + labs(x = "Y", y = "Frequency") + ggtitle(TeX("Original margin of $Y$"))
    expY <- df %>% ggplot(aes(x = Yexp)) + geom_histogram(col = "darkblue", fill = "blue", alpha = 0.3, 
                                                          na.rm = T, show.legend = F) +
      theme_minimal() + labs(x = expression(Y[exp]), y = "Frequency") + ggtitle(TeX("Marginal transformation of $Y$"))
    
    plot_grid(origX, expX, origY, expY, nrow = 2)
  }
  else{
    origjoint <- df %>% ggplot(aes(x = X, y = Y)) + geom_point(na.rm = T) +
      theme_minimal() +
      ggtitle("Original margins")
    expjoint <- df %>% ggplot(aes(x = Xexp, y = Yexp)) + geom_point(na.rm = T) +
      theme_minimal() + labs(x = expression(X[exp]), y = expression(Y[exp])) +
      ggtitle("Standard exponential margins")
    plot_grid(origjoint, expjoint)
  }
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
#' @return An object of S4 class \code{margtransf.class}. This object returns the arguments of the function and extra slot \code{dataexp} containing a matrix with the data on standard margins. To visualise the joint distribution, set \code{joint = T} in the \code{plot} environment. If \code{joint = F} (default), \code{plot} returns the marginal distribution of each variable.
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
#' dataexp <- margtransf(data)
#' 
#' # Plots the marginal distributions of X and Y on original vs standard exponential margins
#' plot(dataexp) 
#' 
#' # Plots the joint distribution of X and Y on original vs standard exponential margins
#' plot(dataexp, joint = T) 
#' 
#' @export
#' 
margtransf <- function(data, qmarg = 0.95){
  if(qmarg < 0 | qmarg > 1){
    stop("Marginal quantile needs to be in [0, 1].")
  }
  result <- margtransf.class(data = data, qmarg = qmarg, dataexp = array())
  nas <- colSums(is.na(data))
  if(any(nas > 0)){
    invisible(sapply(1:length(nas[nas > 0]), function(i){
      warning(paste0("There are ", nas[i], " missing values in margin X", i, ".\n These were removed."))
    }))
  }
  dataunif <- apply(data, 2, empirical_cdf, qmarg = qmarg)
  result@dataexp <- apply(dataunif, 2, qexp)
  return(result)
}  




