#' Goodness of fit of the Angular Dependence function estimates
#' 
#' @name ReturnCurves-package
#' 
#' @description Package
#' 
#' @docType package
#' 
#' @details details \insertCite{MurphyBarltropetal2024;textual}{ReturnCurves}.
#' 
#' @rdname ReturnCurves-package
#' 
#' @aliases ReturnCurves
#' 
#' @references \insertAllCited{}
#' 
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
#' lambda <- adf_est(margdata = margdata, method = "hill")
#' 
#' w_ind <- 31
#' gof <- adf_gof(adf = lambda, w_ind = w_ind)
#' 
#' plot(gof)
#' 
#' \dontrun{
#' # To see the the S4 object's slots
#' str(gof)
#' 
#' # To access the list of vectors
#' gof@@gof
#' }
#' 
#' @export
#' 