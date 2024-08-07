% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_gen_beta.R
\name{gen.data.beta}
\alias{gen.data.beta}
\title{generate data with X and C having beta distributions}
\usage{
gen.data.beta(n, q, B, s2, x.thetas, x.gamma, c.gamma)
}
\arguments{
\item{n}{a positive integer, the sample size}

\item{q}{a number in \link{0,1}, the censoring proportion}

\item{B}{a vector of numbers, parameters in the outcome model}

\item{s2}{a positive number, variance in the outcome model}

\item{x.gamma}{a positive number, gamma parameter for X}

\item{c.gamma}{a positive number, gamma parameter for C}
}
\value{
a list of the following data frames:
\itemize{
\item{\code{datf}: the full data set with Y, X, C, Z}
\item{\code{dat0}: the oracle data Y, W, Delta, Z (with no censoring)}
\item{\code{datcc}: the observed data Y, W, Delta, Z}
\item{\code{dat}: the complete cases from the observed data}
}
}
\description{
generate data with X and C having beta distributions
}
