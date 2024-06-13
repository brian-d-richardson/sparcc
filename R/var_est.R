#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' influence function based variance estimator
#'
#' @param dat data frame, the observed data
#' @param theta a numeric vector, outcome parameters (beta, log(s2))
#' @param args a list of additional arguments
#' @param get.S a function of theta, the estimating function
#' @return.se an optional logical indicator for returning a vector of standard
#' errors as opposed to a covariance matrix, defaults to FALSE
#'
#' @return a covariance matrix or vector of standard errors
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

var.est <- function(dat, theta, args, get.S, return.se = F) {

  # invert empirical mean of outer product of S
  S <- get.S(dat = dat, B = head(theta, -1), ls2 = tail(theta, 1),
             args = args, return.sums = F)
  Sigma <- solve(matrix(rowSums(apply(S, 1, function(s) s %*% t(s))),
                        nrow = length(theta)))

  # return either standard errors or covariance matrix
  if (return.se) {
    return(sqrt(diag(Sigma)))
  } else {
    return(Sigma)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' sandwich variance estimator
#'
#' @inheritParams var.est
#'
#' @return a covariance matrix or vector of standard errors
#'
#' @import numDeriv
#'
#' @export
#'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

var.est.sand <- function(dat, theta, args, get.S, return.se = F) {

  n <- nrow(dat)

  # D: empirical mean of derivative of Psi
  D <- numDeriv::jacobian(
    f = function(x) get.S(dat = dat, B = head(x, -1), ls2 = tail(x, 1),
                          args = args, return.sums = T),
    x = theta,
    method = "simple") / -n
  Dinv <- solve(D)

  # B: empirical mean of outer product of Psi
  S <- get.S(dat = dat, B = head(theta, -1), ls2 = tail(theta, 1),
             args = args, return.sums = F)
  Omega <- matrix(rowMeans(apply(S, 1, function(s) s %*% t(s))),
                  nrow = length(theta))

  Sigma <- Dinv %*% Omega %*% t(Dinv) / n

  # return either standard errors or covariance matrix
  if (return.se) {
    return(sqrt(diag(Sigma)))
  } else {
    return(Sigma)
  }
}


