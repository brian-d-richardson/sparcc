#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' variance estimator
#'
#' @param theta
#' @param get.S
#' @param n
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

var.est <- function(dat, theta, args, get.S, return.se = F) {

  # invert empirical mean of outer product of S
  S <- get.S(dat = dat, B = head(theta, -1), s2 = exp(tail(theta, 1)),
             args = args, return.sums = F)
  Sigma <- solve(matrix(rowMeans(apply(S, 1, function(s) s %*% t(s))),
                        nrow = length(theta))) / nrow(dat)

  # return either standard errors or covariance matrix
  if (return.se) {
    return(sqrt(diag(Sigma)))
  } else {
    return(Sigma)
  }
}



