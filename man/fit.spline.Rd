% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spline_fit.R
\name{fit.spline}
\alias{fit.spline}
\title{Fit spline model with censored outcome data}
\usage{
fit.spline(dat, m.knots, deg)
}
\arguments{
\item{dat}{data frame, the observed data}

\item{m.knots}{a positive integer, the number of knots for the spline model}

\item{deg}{a positive integer, the degree of the spline basis functions}
}
\value{
a list including
\itemize{
\item{\code{alpha}: a matrix, unconstrained spline model parameters}
\item{\code{theta}: a matrix, unconstrained spline model parameters}
\item{\code{dens}: a function, the estimated density}
\item{\code{logf}: a function, the log-likelihood}
\item{\code{glogf}: a function, gradient of the log-likelihood}
\item{\code{deg}: a positive integer, degree of the spline basis functions}
\item{\code{knots}: a numeric vector, knots used in the spline model}
\item{\code{wts}: a numeric vector, weights used in the theta-to-alpha transformation}
}
}
\description{
Fit spline model with censored outcome data
}
