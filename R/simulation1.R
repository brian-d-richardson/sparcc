#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' simulation setup 1
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
sim1 <- function(n, q, B2, s2, x.mean.scale, x.shape, c.shape,
                 specify.x.gamma, specify.c.gamma,
                 mx = 50, mc = 15, my = 3, seed) {

  ## for troublehsooting
  #library(devtools); load_all()
  #n = 4000; q = 0.7; x.mean.scale = 1; s2 = 0.81; B2 = 0.5; x.shape = 1.25; c.shape = 2; mx = 50; mc = 15; my = 3; seed = 1
  #specify.x.gamma = T; specify.c.gamma = T
  set.seed(seed)

  ## define parameters
  B <- c(1, B2, 2)
  x.means = x.mean.scale * c(0.5, 1)

  ## generate data
  dat.list <- gen.data(
    n = n, q = q, B = B, s2 = s2,
    x.means = x.means,
    x.shape = x.shape, c.shape = c.shape)
  datf <- dat.list$datf          # full data
  dat0 <- dat.list$dat0          # oracle data
  dat <- dat.list$dat            # observed data
  datcc <- dat.list$datcc        # complete case data
  zs <- sort(unique(dat$Z))      # unique z values

  ## define nuisance densities

  # estimate distribution of X|Z
  if (specify.x.gamma) {

    # estimate gamma parameters at each level of Z
    x.params.hat <- t(vapply(
      X = 0:1,
      FUN.VALUE = numeric(2),
      FUN = function(z) {
        z.ind <- dat$Z == z
        gammaMLE(yi = dat$W[z.ind],
                 si = dat$Delta[z.ind],
                 scale = F)$estimate
      }))

    # define estimated X|Z density
    eta1 <- function(x, z) {
      dgamma(x = x,
             shape = x.params.hat[z + 1, "shape"],
             rate = x.params.hat[z + 1, "rate"])
    }

  } else {

    # estimate exponential parameter at each level of Z
    x.rates.hat <- vapply(
      X = 0:1,
      FUN.VALUE = 0,
      FUN = function(z) {
        z.ind <- dat$Z == z
        mean(dat$Delta[z.ind]) / mean(dat$W[z.ind])
      })

    # define estimated X|Z density
    eta1 <- function(x, z) {
      dexp(x, rate = x.rates.hat[z + 1])
    }
  }

  # estimate distribution of C|Z
  if (specify.c.gamma) {

    # estimate gamma parameters at each level of Z
    c.params.hat <- t(vapply(
      X = 0:1,
      FUN.VALUE = numeric(2),
      FUN = function(z) {
        z.ind <- dat$Z == z
        gammaMLE(yi = dat$W[z.ind],
                 si = 1 - dat$Delta[z.ind],
                 scale = F)$estimate
      }))

    # define estimated C|Z density
    eta2 <- function(c, z) {
      dgamma(x = c,
             shape = c.params.hat[z + 1, "shape"],
             rate = c.params.hat[z + 1, "rate"])
    }

  } else {

    # estimate exponential parameter at each level of Z
    c.rates.hat <- vapply(
      X = 0:1,
      FUN.VALUE = 0,
      FUN = function(z) {
        z.ind <- dat$Z == z
        mean(1 - dat$Delta[z.ind]) / mean(dat$W[z.ind])
      })

    # define estimated C|Z density
    eta2 <- function(c, z) {
      dexp(x = c, rate = c.rates.hat[z + 1])
    }
  }

  ## create quadrature rules

  # X|Z quadrature
  x.nds <- vapply(
    X = zs,
    FUN.VALUE = numeric(mx),
    FUN = function(z) seq(10^-6, max(datf$X[datf$Z == z]), length = mx))

  x.wts <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mx),
    FUN = function(i) eta1(x.nds[,i], zs[i]) / sum(eta1(x.nds[,i], zs[i])))

  # C|Z quadrature
  c.nds <- vapply(
    X = zs,
    FUN.VALUE = numeric(mc),
    FUN = function(z) seq(10^-6, max(datf$C[datf$Z == z]), length = mc))

  c.wts <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mc),
    FUN = function(i) eta2(c.nds[,i], zs[i]) / sum(eta2(c.nds[,i], zs[i])))

  # Y quadrature
  gq <- gauss.quad(n = my, kind = "hermite")
  y.nds <- gq$nodes
  y.wts <- gq$weights

  ## estimate parameters

  # complete case lm to get starting value
  naive.lm <- lm(Y ~ W + Z, data = datcc)

  # complete case
  Bcc <- get.root(dat = dat, score = get.Scc,
                  start = c(naive.lm$coef, log(var(naive.lm$resid))))

  # oracle
  B0 <- get.root(dat = dat0, score = get.Scc,
                 start = Bcc)

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
           x.mean.scale = x.mean.scale,
           x.shape = x.shape, c.shape = c.shape,
           specify.x.gamma = specify.x.gamma,
           specify.c.gamma = specify.c.gamma,
           mx = mx, mc = mc, my = my, seed = seed,
           Bor = B0, Bcc = Bcc, Bml = Bmle, Bsp = Beff)

  return(ret)
}

# mean function mu(X, Z, B) = E(Y | X, Z)
mu <- function(x, z, B) {
  B[1] + B[2]*x + B[3]*z
}

# gradient of mu w.r.t. B
d.mu <- function(x, z, B) {
  cbind(1, x, z)
}

# Y|X, Z density
fy <- function(y, x, z, B, s2) dnorm(x = y, mean = mu(x, z, B), sd = sqrt(s2))

# full data score vector
SF <- function(y, x, z, B, ls2) {
  cbind((y - mu(x, z, B)) * d.mu(x, z, B),
        (y - mu(x, z, B)) ^ 2 - exp(ls2))
}
