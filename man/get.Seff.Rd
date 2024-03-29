% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_funs.R
\name{get.Seff}
\alias{get.Seff}
\title{semiparametric efficient score function}
\usage{
get.Seff(dat, B, s2, args, return.sums = T)
}
\arguments{
\item{dat}{a data frame including columns Y, W, Delta, Z}

\item{B}{a numeric vector, parameters in the outcome model}

\item{s2}{a positive number, variance in the outcome model}

\item{args}{list of additional arguments}

\item{return.sums}{logical indicator for returning sum of scores as opposed
to individual scores, default is TRUE}
}
\value{
semiparametric efficient score function values
}
\description{
semiparametric efficient score function
}
