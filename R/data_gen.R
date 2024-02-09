#' get censoring probability q given X and C ~ gamma
#'
#' @param
#'
#' @return
#'
#' @export
get.q <- function(x.rate, x.shape, c.rate, c.shape) {

  integrate(f = function(C) {
    C ^ (c.shape - 1) *
      exp(-C * c.rate) *
      pgamma(q = C * x.rate, shape = x.shape, lower.tail = F)
  },
  lower = 0,
  upper = Inf)$value *
    c.rate ^ (c.shape) / gamma(c.shape)
}


#' get rate parameter for C given desired censoring proportion
#'
#' @param
#'
#' @return
#'
#' @export
get.c.rate <- function(q, x.rate, x.shape, c.shape) {

  uniroot(f = function(cr) q - get.q(x.rate = x.rate, x.shape = x.shape,
                                     c.rate = cr, c.shape = c.shape),
          interval = c(1e-3, x.rate),
          extendInt = "downX")$root
}


#' generate data
#'
#' @param
#'
#' @return
#'
#' @export
gen.data <- function(n, q, B, s2, x.rate, x.shape, c.rate, c.shape) {

  X <- rgamma(n, shape = x.shape, rate = x.rate)         # censored covariate
  C <- rgamma(n, shape = c.shape, rate = c.rate)         # censoring time
  W <- ifelse(X <= C, X, C)                              # observed covariate
  Delta <- ifelse(X <= C, 1, 0)                          # uncensored indicator
  Y <- rnorm(n, cbind(1, X) %*% B, sd = sqrt(s2))        # outcome
  datf <- data.frame(Y, X, C)                            # full data
  dat0 <- data.frame(Y, W = X, Delta = 1)                # oracle data
  dat <- data.frame(Y, W, Delta)                         # observed data
  datcc <- dat[Delta == 1,]                              # complete case data

  return(list(datf = datf,
              dat0 = dat0,
              dat = dat,
              datcc = datcc))

}
