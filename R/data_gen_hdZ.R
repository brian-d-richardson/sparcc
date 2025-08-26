#' generate data with X and C having beta distributions with 3D Z
#'
#' @inheritParams get.q.beta
#' @param n a positive integer, the sample size
#' @param q a number in [0,1], the censoring proportion
#' @param B a vector of numbers, parameters in the outcome model
#' @param s2 a positive number, variance in the outcome model
#'
#' @return a list of the following data frames:
#' \itemize{
#' \item{`datf`: the full data set with Y, X, C, Z}
#' \item{`dat0`: the oracle data Y, W, Delta, Z (with no censoring)}
#' \item{`datcc`: the observed data Y, W, Delta, Z}
#' \item{`dat`: the complete cases from the observed data}
#' }
#'
#' @export
gen.data.hdZ <- function(n, q, B, s2, x.thetas, x.gamma, c.gamma) {

  # for troubleshooting
  #n <- 5000; q <- 0.8; B <- c(1, -1, 2); s2 <- 0.16;
  #x.thetas <- seq(0.4, 1.0, length = 8)
  #x.gamma <- 4; c.gamma <- 4

  # theta parameters for beta distribution of C|Z
  c.thetas <- vapply(
    X = x.thetas,
    FUN.VALUE = 0,
    FUN = function(xt)
      get.c.param.beta(q = q, x.theta = xt,
                       x.gamma = x.gamma, c.gamma = c.gamma))

  Z <- sample(0:(length(x.thetas) - 1),                  # uncensored covariate
              size = n, replace = T)
  X <- rbeta(n, 1 + x.gamma - x.thetas[Z + 1],           # censored covariate
                1 + x.gamma + x.thetas[Z + 1])
  C <- rbeta(n, 1 + c.gamma - c.thetas[Z + 1],           # censoring time
                1 + c.gamma + c.thetas[Z + 1])
  W <- ifelse(X <= C, X, C)                              # observed covariate
  Delta <- ifelse(X <= C, 1, 0)                          # uncensored indicator
  Y <- rnorm(n, 0, sd = sqrt(s2)) +                      # outcome
    mu(x = X, z = Z, B = B, xz.interaction = F)
  datf <- data.frame(Y, X, C, Z)                         # full data
  dat0 <- data.frame(Y, W = X, Delta = 1, Z)             # oracle data
  dat <- data.frame(Y, W, Delta, Z)                      # observed data
  datcc <- dat[Delta == 1,]                              # complete case data

  return(list(datf = datf,
              dat0 = dat0,
              dat = dat,
              datcc = datcc))

}
