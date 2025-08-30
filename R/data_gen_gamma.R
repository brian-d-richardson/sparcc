#' get right-censoring probability q given X and C ~ gamma
#'
#' @param x.rate a positive number, rate parameter for gamma distribution of X
#' @param x.shape a positive number, shape parameter for gamma distribution of X
#' @param c.rate a positive number, rate parameter for gamma distribution of C
#' @param c.shape a positive number, shape parameter for gamma distribution of C
#'
#' @return a number in [0,1], the right-censoring proportion q
#'
#' @export
get.q.gamma <- function(x.rate, x.shape, c.rate, c.shape) {

  integrate(f = function(C) {
    C ^ (c.shape - 1) *
      exp(-C * c.rate) *
      pgamma(q = C * x.rate, shape = x.shape, lower.tail = F)
  },
  lower = 0,
  upper = Inf)$value *
    c.rate ^ (c.shape) / gamma(c.shape)
}


#' get rate parameter for C given desired right-censoring proportion
#'
#' @inheritParams get.q.gamma
#' @param q a number in [0,1], the right-censoring proportion
#'
#' @return rate parameter for gamma distribution of C
#'
#' @export
get.c.rate <- function(q, x.rate, x.shape, c.shape) {

  uniroot(f = function(cr) q - get.q.gamma(x.rate = x.rate, x.shape = x.shape,
                                           c.rate = cr, c.shape = c.shape),
          interval = c(1e-3, x.rate),
          extendInt = "downX")$root
}


#' generate data with X and C having gamma distributions
#'
#' @inheritParams get.q.gamma
#' @param n a positive integer, the sample size
#' @param q a number in [0,1], the right-censoring proportion
#' @param B a vector of numbers, parameters in the outcome model
#' @param s2 a positive number, variance in the outcome model
#'
#' @return a list of the following data frames:
#' \itemize{
#' \item{`datf`: the full data set with Y, X, C, Z}
#' \item{`dat0`: the oracle data Y, W, Delta, Z (with no right-censoring)}
#' \item{`datcc`: the observed data Y, W, Delta, Z}
#' \item{`dat`: the complete cases from the observed data}
#' }
#'
#' @export
gen.data.gamma <- function(n, q, B, s2, x.means, x.shape, c.shape) {

  # for troubleshooting
  #n <- 5000; q <- 0.8; B <- c(1, -1, 2); s2 <- 0.81; x.means <- c(0.5, 1);
  #x.shape <- 2; c.shape <- 2;

  x.rates <- x.shape / x.means  # rate parameters for gamma distribution of X|Z
  c.rates <- vapply(            # rate parameters for gamma distribution of C|Z
    X = x.rates,
    FUN.VALUE = 0,
    FUN = function(xr)
      get.c.rate(q = q, x.rate = xr,
                 x.shape = x.shape,
                 c.shape = c.shape))

  Z <- rbinom(n, size = 1, prob = 1/2)                   # uncensored covariate
  X <- rgamma(n, shape = x.shape, rate = x.rates[Z + 1]) # censored covariate
  C <- rgamma(n, shape = c.shape, rate = c.rates[Z + 1]) # censoring variable
  W <- ifelse(X <= C, X, C)                              # observed covariate
  Delta <- ifelse(X <= C, 1, 0)                          # uncensored indicator
  Y <- rnorm(n, cbind(1, X, Z) %*% B, sd = sqrt(s2))     # outcome
  datf <- data.frame(Y, X, C, Z)                         # full data
  dat0 <- data.frame(Y, W = X, Delta = 1, Z)             # oracle data
  dat <- data.frame(Y, W, Delta, Z)                      # observed data
  datcc <- dat[Delta == 1,]                              # complete case data

  return(list(datf = datf,
              dat0 = dat0,
              dat = dat,
              datcc = datcc))

}
