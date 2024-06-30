#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' log-likelihood assuming X ~ Beta
#'
#' @param dat data frame, the observed data
#' @param alpha a numeric vector, beta parameters for X|Z
#'
#' @return a matrix of estimating equation values
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log.fx <- function(dat, alpha) {

  w <- dat$W
  delta <- dat$Delta

  alpha1 <- alpha[1]
  alpha2 <- alpha[2]

  ret <- -log(beta(alpha1, alpha2)) +
    delta * ((alpha1 - 1) * log(w) +
               (alpha2 - 1) * log(1 - w)) +
    (1 - delta) * log(Ibeta(x = w, a = alpha1, b = alpha2, lower = F))

  return(ret)

}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' derivative of log-likelihood assuming X ~ Beta
#'
#' @inheritParams get.Sml
#'
#' @param theta a numeric vector, beta parameters for X|Z
#'
#' @return a matrix of estimating equation values
#'
#' @import numDeriv
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d.log.fx <- function(dat, theta, args = list(), return.sums = T) {

  d <- numDeriv::jacobian(
    func = function(aa) log.fx(dat = dat, alpha = aa),
    x = theta,
    method = "simple"
  )

  if (return.sums) {
    return(colSums(d))
  } else {
    return(d)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' log-likelihood assuming X|Z ~ Beta
#'
#' @param dat data frame, the observed data
#' @param alpha a numeric vector, beta parameters for X|Z
#'
#' @return a matrix of estimating equation values
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log.fxz <- function(dat, alpha) {

  w <- dat$W
  delta <- dat$Delta
  z <- dat$Z

  alpha1 <- alpha[1:2][z+1]
  alpha2 <- alpha[3:4][z+1]

  ret <- -log(beta(alpha1, alpha2)) +
    delta * ((alpha1 - 1) * log(w) +
               (alpha2 - 1) * log(1 - w)) +
    (1 - delta) * log(Ibeta(x = w, a = alpha1, b = alpha2, lower = F))

  return(ret)

}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' derivative of log-likelihood assuming X|Z ~ Beta
#'
#' @inheritParams get.Sml
#'
#' @param theta a numeric vector, beta parameters for X|Z
#'
#' @return a matrix of estimating equation values
#'
#' @import numDeriv
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d.log.fxz <- function(dat, theta, args = list(), return.sums = T) {

  d <- numDeriv::jacobian(
    func = function(aa) log.fxz(dat = dat, alpha = aa),
    x = theta,
    method = "simple"
  )

  if (return.sums) {
    return(colSums(d))
  } else {
    return(d)
  }
}
