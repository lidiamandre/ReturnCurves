% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rc_gof.R
\name{rc_gof.class}
\alias{rc_gof.class}
\title{An S4 class to represent the Goodness-of-Fit of the Return Curve estimates}
\usage{
rc_gof.class(retcurve, blocksize, nboot, alpha, gof)
}
\description{
An S4 class to represent the Goodness-of-Fit of the Return Curve estimates
}
\section{Slots}{

\describe{
\item{\code{retcurve}}{An S4 object of class \code{rc_est.class}.}

\item{\code{blocksize}}{Size of the blocks for the block bootstrap procedure. If \code{1} (default), then a standard bootstrap approach is applied.}

\item{\code{nboot}}{Number of bootstrap samples to be taken. Default is \code{250} samples.}

\item{\code{nangles}}{\loadmathjax{} Number of angles \mjeqn{m}{m} in the interval \mjeqn{(0, \pi/2)}{} \insertCite{MurphyBarltropetal2023}{ReturnCurves}. Default is \code{150} angles.}

\item{\code{alpha}}{Significance level to compute the \mjeqn{(1-\alpha)}{}\% confidence intervals. Default is \code{0.05}.}

\item{\code{gof}}{A list containing the median of the empirical probability and the lower and upper bound of the confidence interval.}
}}

\keyword{internal}
