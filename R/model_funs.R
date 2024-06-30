#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' linear mean function mu(X, Z, B) = E(Y | X, Z)
#'
#' @param x a numeric vector, covariate values
#' @param z a numeric vector, covariate values
#' @param B a numeric vector, parameter value
#' @param xz.interaction logical indicator for including X*Z interaction
#'
#' @return a numeric vector, mean values
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mu <- function(x, z, B, xz.interaction) {
  if (xz.interaction) {
    B[1] + B[2]*x + B[3]*z + B[4]*x*z
  } else {
    B[1] + B[2]*x + B[3]*z
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' gradient of linear mean function mu w.r.t. B
#'
#' @inheritParams mu
#'
#' @return a numeric matrix, the gradient of mu w.r.t. B
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d.mu <- function(x, z, B, xz.interaction) {
  if (xz.interaction) {
    cbind(1, x, z, x*z)
  } else {
    cbind(1, x, z)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' normal density for Y|X,Z
#'
#' @inheritParams mu
#'
#' @param y a numeric vector, outcome values
#' @param s2 a positive number, the variance of Y|X,Z
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fy <- function(y, x, z, B, s2, xz.interaction) {
  dnorm(x = y,
        mean = mu(x, z, B, xz.interaction = xz.interaction),
        sd = sqrt(s2))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Full data score function
#'
#' @inheritParams mu
#'
#' @param ls2 a number, the log of the variance of Y|X,Z
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SF <- function(y, x, z, B, ls2, xz.interaction) {
  cbind((y - mu(x, z, B, xz.interaction = xz.interaction)) *
          d.mu(x, z, B, xz.interaction = xz.interaction) / exp(ls2),
        0.5 * ( (y - mu(x, z, B, xz.interaction = xz.interaction)) ^ 2 /
                  exp(ls2) - 1 ))
}


