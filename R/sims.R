#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' simulation setup 0
#'
#' @inheritParams get.q
#' @inheritParams gen.data
#'
#' @param B2 a number, the coefficient for `X` in the outcome model
#' @param s2 a number, the variance in the outcome model
#' @param specify.x.gamma logical, an indicator of whether X should be estimated
#' as gamma (`T`) or exponential (`F`)
#' @param specify.c.gamma logical, an indicator of whether C should be estimated
#' as gamma (`T`) or exponential (`F`)
#' @param mx a positive number, the number of nodes in quadrature grid for X
#' @param mc a positive number, the number of nodes in quadrature grid for C
#' @param my a positive number, the number of nodes in quadrature grid for Y
#' @param seed a positive integer, the seed value for random number generation
#'
#' @importFrom statmod gauss.quad
#' @importFrom conf gammaMLE
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sim0 <- function(n, q, B2, s2, x.mean, x.shape, c.shape,
                 specify.x.gamma, specify.c.gamma,
                 mx = 100, mc = 15, my = 3, seed) {

  ## for troublehsooting
  #library(devtools); load_all()
  #n = 10000; q = 0.8; x.mean = 1; s2 = 3; x.shape = 1.2; c.shape = 2; mx = 100; mc = 15; my = 3; seed = 1
  #specify.x.gamma = T; specify.c.gamma = T
  set.seed(seed)

  ## define parameters
  B <- c(1, B2)              # outcome model parameters
  x.rate <- x.shape / x.mean # rate parameter for gamma distribution of X
  c.rate <- get.c.rate(      # rate parameter for gamma distribution of C
    q = q,
    x.rate = x.rate,
    x.shape = x.shape,
    c.shape = c.shape)

  ## generate data
  dat.list <- gen.data(n = n, q = q, B = B, s2 = s2,
                       x.rate = x.rate, x.shape = x.shape,
                       c.rate = c.rate, c.shape = c.shape)
  datf <- dat.list$datf          # full data
  dat0 <- dat.list$dat0          # oracle data
  dat <- dat.list$dat            # observed data
  datcc <- dat.list$datcc        # complete case data

  ## define nuisance densities

  # distribution of X
  if (specify.x.gamma) {
    eta1 <- function(x) dgamma(x = x, shape = x.shape, rate = x.rate)
  } else {
    eta1 <- function(x) dexp(x, rate = x.mean)
  }

  # distribution of C
  if (specify.c.gamma) {
    eta2 <- function(x) dgamma(x = x, shape = c.shape, rate = c.rate)
  } else {
    eta2 <- function(x) dexp(x, rate = c.rate / c.shape)
  }

  x.upper <- max(datf$X)  # oracle max
  c.upper <- max(datf$C)

  ## create quadrature rules

  # X quadrature
  x.nds <- seq(10^-5, x.upper, length = mx)
  x.wts <- eta1(x.nds) / sum(eta1(x.nds))

  # C quadrature
  c.nds <- seq(10^-5, c.upper, length = mc)
  c.wts <- eta2(c.nds) / sum(eta2(c.nds))

  # Y quadrature
  gq <- statmod::gauss.quad(n = my, kind = "hermite")
  y.nds <- gq$nodes
  y.wts <- gq$weights

  ## estimate parameters

  # complete case lm to get starting value
  naive.lm <- lm(Y ~ W, data = datcc)

  # oracle
  B0 <- get.root(dat = dat0, score = get.Scc,
                 start = c(naive.lm$coef, log(var(naive.lm$resid))))

  # complete case
  Bcc <- get.root(dat = dat, score = get.Scc,
                  start = c(naive.lm$coef, log(var(naive.lm$resid))))

  # MLE
  Bmle <- get.root(dat = dat, score = get.Sml, start = Bcc,
                   args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                               x.nds = x.nds, x.wts = x.wts))

  # semiparametric efficient score
  Beff <- get.root(dat = dat, score = get.Seff, start = Bcc,
                   args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                               eta1 = eta1,
                               x.nds = x.nds, x.wts = x.wts,
                               c.nds = c.nds, c.wts = c.wts,
                               y.nds = y.nds, y.wts = y.wts))

  # return setup parameters and estimates
  ret <- c(n = n, q = q, B2 = B2, s2 = s2,
           x.mean = x.mean, x.shape = x.shape, c.shape = c.shape,
           specify.x.gamma = specify.x.gamma,
           specify.c.gamma = specify.c.gamma,
           mx = mx, mc = mc, my = my, seed = seed,
           Bor = B0, Bcc = Bcc, Bml = Bmle, Bsp = Beff)

  return(ret)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' simulation setup 1
#'
#' @inheritParams sim0
#' @inheritParams gen.data
#'
#' @importFrom statmod gauss.quad
#' @importFrom conf gammaMLE
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sim1 <- function(n, q, B2, s2, x.mean, x.shape, c.shape,
                 specify.x.gamma, specify.c.gamma,
                 mx = 100, mc = 15, my = 3, seed) {

  ## for troublehsooting
  #library(devtools); load_all()
  #n = 10000; q = 0.8; x.mean = 0.5; s2 = 1.1; x.shape = 1.2; c.shape = 1.2; mx = 100; mc = 15; my = 2; seed = 1

  set.seed(seed)

  ## define parameters
  B <- c(1, B2)              # outcome model parameters
  x.rate <- x.shape / x.mean # rate parameter for gamma distribution of X
  c.rate <- get.c.rate(      # rate parameter for gamma distribution of C
    q = q,
    x.rate = x.rate,
    x.shape = x.shape,
    c.shape = c.shape)

  ## generate data
  dat.list <- gen.data(n = n, q = q, B = B, s2 = s2,
                       x.rate = x.rate, x.shape = x.shape,
                       c.rate = c.rate, c.shape = c.shape)
  datf <- dat.list$datf          # full data
  dat0 <- dat.list$dat0          # oracle data
  dat <- dat.list$dat            # observed data
  datcc <- dat.list$datcc        # complete case data

  ## estimate nuisance densities

  # estimate distribution of X
  if (specify.x.gamma) {
    x.param.hat <- gammaMLE(yi = dat$W, si = dat$Delta, scale = F)$estimate
    eta1 <- function(x)
      dgamma(x = x, shape = x.param.hat["shape"], rate = x.param.hat["rate"])
  } else {
    x.rate.hat <- mean(dat$Delta) / mean(dat$W)
    eta1 <- function(x) dexp(x, rate = x.rate.hat)
  }

  # estimate distribution of C
  if (specify.c.gamma) {
    c.param.hat <- gammaMLE(yi = dat$W, si = 1 - dat$Delta, scale = F)$estimate
    eta2 <- function(x)
      dgamma(x = x, shape = c.param.hat["shape"], rate = c.param.hat["rate"])
  } else {
    c.rate.hat <- mean(1 - dat$Delta) / mean(dat$W)
    eta2 <- function(x) dexp(x, rate = c.rate.hat)
  }

  x.upper <- max(datf$X)  # oracle max
  c.upper <- max(datf$C)
  #x.upper <- sum(1/(1:n)) / x.rate.hat # estimated expected max
  #c.upper <- sum(1/(1:n)) / c.rate.hat
  #x.upper <- qexp((n-1)/n, rate = x.rate.hat)  # estimated (n-1)/n quantile
  #c.upper <- qexp((n-1)/n, rate = c.rate.hat)

  ## create quadrature rules

  # X quadrature
  x.nds <- seq(10^-5, x.upper, length = mx)
  x.wts <- eta1(x.nds) / sum(eta1(x.nds))

  # C quadrature
  c.nds <- seq(10^-5, c.upper, length = mc)
  c.wts <- eta2(c.nds) / sum(eta2(c.nds))

  # Y quadrature
  gq <- statmod::gauss.quad(n = my, kind = "hermite")
  y.nds <- gq$nodes
  y.wts <- gq$weights

  ## estimate parameters

  # complete case lm to get starting value
  naive.lm <- lm(Y ~ W, data = datcc)

  # oracle
  B0 <- get.root(dat = dat0, score = get.Scc,
                 start = c(naive.lm$coef, log(var(naive.lm$resid))))

  # complete case
  Bcc <- get.root(dat = dat, score = get.Scc,
                  start = c(naive.lm$coef, log(var(naive.lm$resid))))

  # MLE
  Bmle <- get.root(dat = dat, score = get.Sml, start = Bcc,
                   args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                               x.nds = x.nds, x.wts = x.wts))

  # semiparametric efficient score
  Beff <- get.root(dat = dat, score = get.Seff, start = Bcc,
                   args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                               eta1 = eta1,
                               x.nds = x.nds, x.wts = x.wts,
                               c.nds = c.nds, c.wts = c.wts,
                               y.nds = y.nds, y.wts = y.wts))

  # return setup parameters and estimates
  ret <- c(n = n, q = q, B2 = B2, s2 = s2,
           x.mean = x.mean, x.shape = x.shape, c.shape = c.shape,
           specify.x.gamma = specify.x.gamma,
           specify.c.gamma = specify.c.gamma,
           mx = mx, mc = mc, my = my, seed = seed,
           Bor = B0, Bcc = Bcc, Bml = Bmle, Bsp = Beff)

  return(ret)

}


# mean function mu(X, B) = E(Y | X)
mu <- function(x, B) {
  B[1] + B[2]*x
}

# gradient of mu w.r.t. B
d.mu <- function(x, B) {
  cbind(1, x)
}

# Y density
fy <- function(y, x, B, s2) dnorm(x = y, mean = mu(x, B), sd = sqrt(s2))

# full data score vector
SF <- function(y, x, B, s2) {
  cbind((y - mu(x, B)) * d.mu(x, B),
        (y - mu(x, B)) ^ 2 - s2)
}

