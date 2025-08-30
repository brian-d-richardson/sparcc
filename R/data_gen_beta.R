#' get right-censoring probability q given X and C ~ Beta
#'
#' @param x.theta a positive number, theta parameter for X
#' @param x.gamma a positive number, gamma parameter for X
#' @param c.theta a positive number, theta parameter for C
#' @param c.gamma a positive number, gamma parameter for C
#'
#' @importFrom zipfR Ibeta
#'
#' @return a number in [0,1], the right-censoring proportion q
#'
#' @export
get.q.beta <- function(x.theta, x.gamma, c.theta, c.gamma) {

  # for troubleshooting
  #x.theta <- 1; x.gamma <- 4; c.theta <- 1; c.gamma <- 4

  integrate(
    f = function(x) dbeta(x, 1 + x.gamma - x.theta, 1 + x.gamma + x.theta) *
      zipfR::Ibeta(x = x,
                   lower = T,
                   a = 1 + c.gamma - c.theta,
                   b = 1 + c.gamma + c.theta,
                   log = F),
    lower = 0,
    upper = 1)[[1]] /
     beta(1 + c.gamma - c.theta, 1 + c.gamma + c.theta)
}


#' get beta parameters for C given desired right-censoring proportion
#'
#' @inheritParams get.q.beta
#' @param q a number in [0,1], the right-censoring proportion
#'
#' @return a positive number, the theta parameter for C
#'
#' @export
get.c.param.beta <- function(q, x.theta, x.gamma, c.gamma) {

  uniroot(f = function(ct) q - get.q.beta(x.theta = x.theta,
                                          x.gamma = x.gamma,
                                          c.theta = ct,
                                          c.gamma = c.gamma),
          interval = c(-c.gamma, c.gamma),
          extendInt = "downX")$root
}


#' generate data with X and C having beta distributions
#'
#' @inheritParams get.q.beta
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
gen.data.beta <- function(n, q, B, s2, x.thetas, x.gamma, c.gamma) {

  # for troubleshooting
  #n <- 5000; q <- 0.8; B <- c(1, -1, 2); s2 <- 0.16; x.thetas <- c(0.5, 1);
  #x.gamma <- 4; c.gamma <- 4

  # theta parameters for beta distribution of C|Z
  c.thetas <- vapply(
    X = x.thetas,
    FUN.VALUE = 0,
    FUN = function(xt)
      get.c.param.beta(q = q, x.theta = xt,
                       x.gamma = x.gamma, c.gamma = c.gamma))

  Z <- rbinom(n, size = 1, prob = 1/2)                   # uncensored covariate
  X <- rbeta(n, 1 + x.gamma - x.thetas[Z + 1],           # censored covariate
                1 + x.gamma + x.thetas[Z + 1])
  C <- rbeta(n, 1 + c.gamma - c.thetas[Z + 1],           # censoring value
                1 + c.gamma + c.thetas[Z + 1])
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
