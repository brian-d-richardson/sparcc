#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' simulation setup 4
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
#' @import DALSM
#' @import dplyr
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sim5 <- function(n, q, seed) {

  ## for troublehsooting
  #rm(list = ls()); library(devtools); load_all()
  #n = 8000; q = 0.8; seed <- 1

  # define parameters -------------------------------------------------------

  B <- c(1, 10, 2)
  s2 <- 1
  x.thetas <- 0.5 * c(-1, 1)
  x.gamma <- 1
  c.gamma <- 2
  mx <- 40
  mc <- 40
  my <- 5

  # generate data -----------------------------------------------------------

  set.seed(seed)
  dat.list <- gen.data.beta(
    n = n, q = q, B = B, s2 = s2,
    x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma)

  datf <- dat.list$datf          # full data
  dat0 <- dat.list$dat0          # oracle data
  dat <- dat.list$dat            # observed data
  datcc <- dat.list$datcc        # complete case data
  zs <- sort(unique(dat$Z))      # unique z values

  # estimate nuisance distributions -----------------------------------------

  ## X|Z quadrature nodes
  x.nds <- vapply(
    X = zs,
    FUN.VALUE = numeric(mx),
    FUN = function(z) seq(1E-6, 1-1E-6, length = mx))

  ## X|Z correctly modeled as conditional beta
  x.params.hat.correct <- t(vapply(
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

  eta1.correct <- function(x, z) {
    dbeta(x = x,
          shape1 = x.params.hat.correct[z + 1, "shape1"],
          shape2 = x.params.hat.correct[z + 1, "shape2"])
  }

  x.wts.correct <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mx),
    FUN = function(i) eta1.correct(x.nds[,i], zs[i]) /
      sum(eta1.correct(x.nds[,i], zs[i])))

  ## X|Z incorrectly modeled as marginal beta
  x.fit.wrong <- dat %>%
    mutate(left = W,
           right = ifelse(Delta == 1, W, NA)) %>%
    dplyr::select(left, right) %>%
    fitdistrplus::fitdistcens(distr = "beta")
  x.params.hat.wrong <- x.fit.wrong$estimate

  eta1.wrong <- function(x, z) {
    dbeta(x = x,
          shape1 = x.params.hat.wrong["shape1"],
          shape2 = x.params.hat.wrong["shape2"])
  }

  x.wts.wrong <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mx),
    FUN = function(i) eta1.wrong(x.nds[,i], zs[i]) /
      sum(eta1.wrong(x.nds[,i], zs[i])))

  ## nonparametric estimated distribution of X|Z
  x.dens.hat.nonpar <- lapply(
    as.list(zs),
    function(z) densityLPS(Dens1d(
      y = dat$W[dat$Z == z],
      event = dat$Delta[dat$Z == z],
      ymin = 0, ymax = 1))
  )

  eta1.nonpar <- function(x, z) {
    eta <- numeric(length(x))
    for (zz in unique(z)) {
      eta[z == zz] <- x.dens.hat.nonpar[[zz + 1]]$ddist(x[z == zz])
    }
    return(eta)
  }

  x.wts.nonpar <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mx),
    FUN = function(i) eta1.nonpar(x.nds[,i], zs[i]) /
      sum(eta1.nonpar(x.nds[,i], zs[i])))

  ## C|Z quadrature nodes
  c.nds <- vapply(
    X = zs,
    FUN.VALUE = numeric(mc),
    FUN = function(z) seq(1E-6, 1-1E-6, length = mc))

  ## C|Z correctly modeled as conditional beta
  c.params.hat.correct <- t(vapply(
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

  eta2.correct <- function(c, z) {
    dbeta(x = c,
          shape1 = c.params.hat.correct[z + 1, "shape1"],
          shape2 = c.params.hat.correct[z + 1, "shape2"])
  }

  c.wts.correct <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mc),
    FUN = function(i) eta2.correct(c.nds[,i], zs[i]) /
      sum(eta2.correct(c.nds[,i], zs[i])))

  ## C|Z incorrectly modeled as marginal beta
  c.fit.wrong <- dat %>%
    mutate(left = W,
           right = ifelse(Delta == 0, W, NA)) %>%
    dplyr::select(left, right) %>%
    fitdistrplus::fitdistcens(distr = "beta")
  c.params.hat.wrong <- c.fit.wrong$estimate

  eta2.wrong <- function(c, z) {
    dbeta(x = c,
          shape1 = c.params.hat.wrong["shape1"],
          shape2 = c.params.hat.wrong["shape2"])
  }

  c.wts.wrong <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mc),
    FUN = function(i) eta2.wrong(c.nds[,i], zs[i]) /
      sum(eta2.wrong(c.nds[,i], zs[i])))

  ## nonparametric estimated distribution of C|Z
  c.dens.hat.nonpar <- lapply(
    as.list(zs),
    function(z) densityLPS(Dens1d(
      y = dat$W[dat$Z == z],
      event = 1 - dat$Delta[dat$Z == z],
      ymin = 0, ymax = 1))
  )

  eta2.nonpar <- function(c, z) {
    eta <- numeric(length(c))
    for (zz in unique(z)) {
      eta[z == zz] <- c.dens.hat.nonpar[[zz + 1]]$ddist(c[z == zz])
    }
    return(eta)
  }

  c.wts.nonpar <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mc),
    FUN = function(i) eta2.nonpar(c.nds[,i], zs[i]) /
      sum(eta2.nonpar(c.nds[,i], zs[i])))

  ## Y|X,Z quadrature nodes
  gq <- gauss.quad(n = my, kind = "hermite")
  y.nds <- gq$nodes
  y.wts <- gq$weights

  # estimate beta -----------------------------------------------------------

  ## complete case linear model to get starting value
  naive.lm <- lm(Y ~ W + Z, data = datcc)

  ## complete case
  B.cc <- get.root(dat = dat, score = get.Scc,
                   start = c(naive.lm$coef, log(var(naive.lm$resid))))
  V.cc <- var.est(dat = datcc, theta = B.cc, args = list(),
                  get.S = get.Scc, return.se = T)

  ## oracle
  B.or <- get.root(dat = dat0, score = get.Scc, start = B.cc)
  V.or <- var.est(dat = dat0, theta = B.or, args = list(),
                  get.S = get.Scc, return.se = T)

  ## MLE (X|Z correct)
  mle.args <- list(mu = mu, d.mu = d.mu, SF = SF, fy = fy, x.nds = x.nds)
  B.ml.1 <- get.root(dat = dat, score = get.Sml, start = B.cc,
                     args = append(mle.args, list("x.wts" = x.wts.correct)))
  V.ml.1 <- var.est(dat = dat, theta = B.ml.1,
                    args = append(mle.args, list("x.wts" = x.wts.correct)),
                    get.S = get.Sml, return.se = T)

  ## MLE (X|Z incorrect)
  B.ml.0 <- get.root(dat = dat, score = get.Sml, start = B.cc,
                     args = append(mle.args, list("x.wts" = x.wts.wrong)))
  V.ml.0 <- var.est(dat = dat, theta = B.ml.0,
                    args = append(mle.args, list("x.wts" = x.wts.wrong)),
                    get.S = get.Sml, return.se = T)

  ## MLE (X|Z nonparametric)
  B.ml.2 <- get.root(dat = dat, score = get.Sml, start = B.cc,
                     args = append(mle.args, list("x.wts" = x.wts.nonpar)))
  V.ml.2 <- var.est(dat = dat, theta = B.ml.2,
                    args = append(mle.args, list("x.wts" = x.wts.nonpar)),
                    get.S = get.Sml, return.se = T)

  ## semiparametric (X|Z, C|Z correct)
  sp.args <- list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                  x.nds = x.nds, c.nds = c.nds,
                  y.nds = y.nds, y.wts = y.wts)
  B.sp.11 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.correct,
                                                  "c.wts" = c.wts.correct,
                                                  "eta1" = eta1.correct)))
  V.sp.11 <- var.est(dat = dat, theta = B.sp.11,
                     args = append(sp.args, list("x.wts" = x.wts.correct,
                                              "c.wts" = c.wts.correct,
                                              "eta1" = eta1.correct)),
                     get.S = get.Seff, return.se = T)

  ## semiparametric (X|Z correct, C|Z incorrect)
  B.sp.10 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.correct,
                                                  "c.wts" = c.wts.wrong,
                                                  "eta1" = eta1.correct)))
  V.sp.10 <- var.est(dat = dat, theta = B.sp.10,
                     args = append(sp.args, list("x.wts" = x.wts.correct,
                                                 "c.wts" = c.wts.wrong,
                                                 "eta1" = eta1.correct)),
                     get.S = get.Seff, return.se = T)

  ## semiparametric (X|Z incorrect, C|Z correct)
  B.sp.01 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.wrong,
                                                  "c.wts" = c.wts.correct,
                                                  "eta1" = eta1.wrong)))
  V.sp.01 <- var.est(dat = dat, theta = B.sp.01,
                     args = append(sp.args, list("x.wts" = x.wts.wrong,
                                                 "c.wts" = c.wts.correct,
                                                 "eta1" = eta1.wrong)),
                     get.S = get.Seff, return.se = T)

  ## semiparametric (X|Z, C|Z incorrect)
  B.sp.00 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.wrong,
                                                  "c.wts" = c.wts.wrong,
                                                  "eta1" = eta1.wrong)))
  V.sp.00 <- var.est(dat = dat, theta = B.sp.00,
                     args = append(sp.args, list("x.wts" = x.wts.wrong,
                                                 "c.wts" = c.wts.wrong,
                                                 "eta1" = eta1.wrong)),
                     get.S = get.Seff, return.se = T)

  ## semiparametric (X|Z, C|Z nonparametric)
  B.sp.22 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.nonpar,
                                                  "c.wts" = c.wts.nonpar,
                                                  "eta1" = eta1.nonpar)))
  V.sp.22 <- var.est(dat = dat, theta = B.sp.22,
                     args = append(sp.args, list("x.wts" = x.wts.nonpar,
                                                 "c.wts" = c.wts.nonpar,
                                                 "eta1" = eta1.nonpar)),
                     get.S = get.Seff, return.se = T)

  # return setup parameters, estimates, and standard errors (length 83)
  ret <- c(n = n, q = q, seed = seed,
           B.or = B.or, B.cc = B.cc,
           B.ml.0 = B.ml.0, B.ml.1 = B.ml.1, B.ml.2 = B.ml.2,
           B.sp.00 = B.sp.00, B.sp.01 = B.sp.01, B.sp.10 = B.sp.10,
           B.sp.11 = B.sp.11, B.sp.22 = B.sp.22,
           V.or = V.or, V.cc = V.cc,
           V.ml.0 = V.ml.0, V.ml.1 = V.ml.1, V.ml.2 = V.ml.2,
           V.sp.00 = V.sp.00, V.sp.01 = V.sp.01, V.sp.10 = V.sp.10,
           V.sp.11 = V.sp.11, V.sp.22 = V.sp.22)

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




