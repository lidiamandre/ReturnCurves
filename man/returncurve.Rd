% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/returncurve.R
\docType{methods}
\name{rc_est}
\alias{rc_est}
\title{Return Curve Estimation}
\usage{
rc_est(
  data,
  w,
  p,
  method = c("hill", "cl"),
  q = 0.95,
  k = 7,
  constrained = "no"
)
}
\arguments{
\item{data}{matrix that contains the data, in standard exponential margins}

\item{w}{sequence of angles between 0 and 1; default set to a vector of 101 equally spaced angles}

\item{p}{probability for the return curve}

\item{method}{method to be used in the estimation of the angular dependence function: "hill" to use the Hill estimator, "cl" for the composite likelihood estimator}

\item{q}{quantile to be used for the Hill estimator and/or the Heffernan and Tawn conditional extremes model; default set to 0.95}

\item{k}{polynomial degree for the Bernstein-Bezier polynomials used in the estimation of the angular dependence function using the composite likelihood method; default set to 7}

\item{constrained}{indicates whether or not to incorporate knowledge of the conditional extremes parameters; default set to "no"}
}
\value{
return curve estimation
}
\description{
computes the return curve estimation... to do
}
\details{
\loadmathjax This function estimates the return curve given by \mjdeqn{RC(p):=\brace(x, y) \in \mathbb{R}^2: \text{Pr}(X>x, Y>y)=p\brace.} ... talk about how it connects to the estimation of the adf and the methods used, reference the est_lamb function
}
\examples{
library(ReturnCurves)

}
\references{
to do
}