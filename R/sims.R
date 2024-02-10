#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' simulation setup 1
#'
#' @inheritParams get.q
#' @inheritParams gen.data
#'
#' @param mx a positive number, the number of nodes in quadrature grid for X
#' @param mc a positive number, the number of nodes in quadrature grid for C
#' @param my a positive number, the number of nodes in quadrature grid for Y
#' @param seed a positive integer, the seed value for random number generation
#'
#' @importFrom statmod gauss.quad
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sim1 <- function(n, q, x.shape = 1, c.shape = 1,
                 mx = 100, mc = 15, my = 2, seed) {

  ## for troublehsooting
  #library(devtools); load_all()
  #n = 10000; q = 0.8; x.shape = 1; c.shape = 1; mx = 100; mc = 15; my = 2; seed = 1

  set.seed(seed)

  ## define parameters
  B <- c(1, 2)               # outcome model parameters
  s2 <- 1.1                  # Var(Y|X,Z)
  x.mean <- 0.25
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

  ## define densities

  # X density
  eta1 <- function(x) dexp(x, rate = x.rate / x.shape)

  # C density
  eta2 <- function(c) dexp(c, rate = c.rate / c.shape)

  ## create quadrature rules

  # X quadrature
  x.nds <- seq(10^-6, max(datf$X), length = mx)
  x.wts <- eta1(x.nds) / sum(eta1(x.nds))

  # C quadrature
  c.nds <- seq(10^-6, max(datf$C), length = mc)
  c.wts <- eta2(c.nds) / sum(eta2(c.nds))

  # Y quadrature
  gq <- statmod::gauss.quad(n = my, kind = "hermite")
  y.nds <- gq$nodes
  y.wts <- gq$weights

  ## estimate parameters

  # oracle
  B0 <- get.root(dat = dat0, score = get.Scc, start = c(0, 0, 0))

  # complete case
  Bcc <- get.root(dat = dat, score = get.Scc, start = c(0, 0, 0))

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
  ret <- c(n = n, q = q, x.shape = x.shape, c.shape = c.shape,
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

