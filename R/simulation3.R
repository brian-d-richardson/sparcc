#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' simulation setup 3
#'
#' @inheritParams get.q
#' @inheritParams gen.data
#'
#' @param B2 a number, the coefficient for `X` in the outcome model
#' @param s2 a number, the variance in the outcome model
#' @param mx a positive number, the number of nodes in quadrature grid for X
#' @param mc a positive number, the number of nodes in quadrature grid for C
#' @param my a positive number, the number of nodes in quadrature grid for Y
#' @param seed a positive integer, the seed value for random number generation
#'
#' @importFrom statmod gauss.quad
#' @importFrom fitdistrplus fitdistcens
#' @import dplyr
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sim3 <- function(n, q, B2, s2, x.theta.scale, x.gamma, c.gamma,
                 mx = 40, mc = 40, my = 4, seed) {

  ## for troublehsooting
  #library(devtools); load_all()
  #n = 8000; q = 0.8; B2 <- 10; s2 <- 1; x.theta.scale = 0.5; x.gamma <- 2; c.gamma <- 2;
  #mx <- 40; mc <- 40; my <- 5;

  set.seed(seed)

  ## define parameters
  B <- c(1, B2, 2)
  x.thetas <- x.theta.scale * c(-1, 1)

  dat.list <- gen.data.beta(
    n = n, q = q, B = B, s2 = s2,
    x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma)

  datf <- dat.list$datf          # full data
  dat0 <- dat.list$dat0          # oracle data
  dat <- dat.list$dat            # observed data
  datcc <- dat.list$datcc        # complete case data
  zs <- sort(unique(dat$Z))      # unique z values

  ## estimate nuisance distributions

  # estimate beta parameters at each level of Z
  x.params.hat <- t(vapply(
    X = 0:1,
    FUN.VALUE = numeric(2),
    FUN = function(z) {
      est <- dat %>%
        filter(Z == z) %>%
        mutate(left = W,
               right = ifelse(Delta == 1, W, NA)) %>%
        dplyr::select(left, right) %>%
        fitdistrplus::fitdistcens(distr = "beta")
      return(est$estimate)
    }))

  # define estimated X|Z density
  eta1 <- function(x, z) {
    dbeta(x = x,
          shape1 = x.params.hat[z + 1, "shape1"],
          shape2 = x.params.hat[z + 1, "shape2"])
  }


  # estimate beta parameters at each level of Z
  c.params.hat <- t(vapply(
    X = 0:1,
    FUN.VALUE = numeric(2),
    FUN = function(z) {
      est <- dat %>%
        filter(Z == z) %>%
        mutate(left = W,
               right = ifelse(Delta == 0, W, NA)) %>%
        dplyr::select(left, right) %>%
        fitdistrplus::fitdistcens(distr = "beta")
      return(est$estimate)
    }))

  # define estimated C|Z density
  eta2 <- function(c, z) {
    dbeta(x = c,
          shape1 = c.params.hat[z + 1, "shape1"],
          shape2 = c.params.hat[z + 1, "shape2"])
  }

  ## create quadrature rules

  # X|Z quadrature
  x.nds <- vapply(
    X = zs,
    FUN.VALUE = numeric(mx),
    FUN = function(z) seq(1E-6, 1-1E-6, length = mx))

  x.wts <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mx),
    FUN = function(i) eta1(x.nds[,i], zs[i]) / sum(eta1(x.nds[,i], zs[i])))

  # C|Z quadrature
  c.nds <- vapply(
    X = zs,
    FUN.VALUE = numeric(mc),
    FUN = function(z) seq(1E-6, 1-1E-6, length = mc))

  c.wts <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mc),
    FUN = function(i) eta2(c.nds[,i], zs[i]) / sum(eta2(c.nds[,i], zs[i])))

  # Y quadrature
  gq <- gauss.quad(n = my, kind = "hermite")
  y.nds <- gq$nodes
  y.wts <- gq$weights

  ## estimate beta

  # complete case lm to get starting value
  naive.lm <- lm(Y ~ W + Z, data = datcc)

  # complete case
  Bcc <- get.root(dat = dat, score = get.Scc,
                  start = c(naive.lm$coef, log(var(naive.lm$resid))))

  # oracle
  B0 <- get.root(dat = dat0, score = get.Scc, start = Bcc)

  # MLE
  mle.args <- list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                   x.nds = x.nds, x.wts = x.wts)
  Bmle <- get.root(dat = dat, score = get.Sml, start = Bcc,
                   args = mle.args)

  # semiparametric efficient score
  sp.args <- list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                  eta1 = eta1,
                  x.nds = x.nds, x.wts = x.wts,
                  c.nds = c.nds, c.wts = c.wts,
                  y.nds = y.nds, y.wts = y.wts)
  Beff <- get.root(dat = dat, score = get.Seff, start = Bcc,
                   args = sp.args)

  # return setup parameters and estimates
  ret <- c(n = n, q = q, B2 = B2, s2 = s2,
           x.theta.scale = x.theta.scale,
           x.gamma = x.gamma, c.gamma = c.gamma,
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
SF <- function(y, x, z, B, s2) {
  cbind((y - mu(x, z, B)) * d.mu(x, z, B),
        (y - mu(x, z, B)) ^ 2 - s2)
}




