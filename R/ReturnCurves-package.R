#' Estimation of Return Curves
#' 
#' @description \loadmathjax{}
#' Implements the estimation of the \mjeqn{p}{p}-probability return curve \insertCite{MurphyBarltropetal2023}{ReturnCurves}, 
#' as well as a pointwise and smooth estimation of the angular dependence function \insertCite{WadsworthTawn2013}{ReturnCurves}.
#' 
#' @section Available functions: 
#' \code{\link{adf_est}}: Estimation of the Angular Dependence Function (ADF)
#' 
#' \code{\link{adf_gof}}: Goodness of fit of the Angular Dependence Function estimates
#' 
#' \code{\link{margtransf}}: Marginal Transformation
#' 
#' \code{\link{rc_est}}: Return Curve estimation
#'
#' \code{\link{rc_gof}}: Goodness of fit of the Return Curve estimates
#' 
#' \code{\link{rc_unc}}: Uncertainty of the Return Curve estimates
#' 
#' @references \insertAllCited{}
#' 
#' @examples
#' library(ReturnCurves)
#'
#' # Generating data for illustration purposes
#' set.seed(321)
#' data <- cbind(rnorm(1000), rnorm(1000))
#' 
#' n <- dim(data)[1]
#' 
#' # Marginal Transformation
#' margdata <- margtransf(data)
#' 
#' head(margdata@@dataexp)
#' 
#' # Return Curves estimation
#' 
#' prob <- 1/n
#' 
#' retcurve <- rc_est(margdata = margdata, p = prob, method = "hill")
#' 
#' head(retcurve@@rc)
#' 
#' # ADF estimation
#' lambda <- adf_est(margdata = margdata, method = "hill")
#' 
#' head(lambda@@adf)
#' 
#' @export
"_PACKAGE"