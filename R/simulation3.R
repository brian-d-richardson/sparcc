#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' simulation 1: assess various methods with 3D Z
#'
#' @inheritParams get.q.beta
#' @inheritParams gen.data.beta
#'
#' @param ndistinct.Z a positive integer, the number of distinct Z values
#' @param seed a positive integer, the seed value for random number generation
#'
#' @importFrom statmod gauss.quad
#' @importFrom fitdistrplus fitdistcens
#' @importFrom zipfR Ibeta
#' @import dplyr
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sim3 <- function(n, q, ndistinct.Z, seed) {

  ## for troublehsooting
  #rm(list = ls()); library(devtools); load_all(); library(numDeriv)
  #n = 8000; q = 0.8; ndistinct.Z = 4; seed <- 1

  # define parameters -------------------------------------------------------

  # outcome model parameters
  B <- c(1, 10, 2)
  s2 <- 1

  # nuisance model parameters
  x.thetas <- seq(-0.5, 0.5, length = ndistinct.Z)
  x.gamma <- 1
  c.gamma <- 2

  # quadrature rule parameters
  mx <- 40
  mc <- 40
  my <- 5

  # spline estimation parameters
  m.knots <- 5
  deg <- 3

  # generate data -----------------------------------------------------------

  set.seed(seed)
  dat.list <- gen.data.hdZ(
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
  st <- Sys.time()
  x.params.hat.correct <- t(vapply(
    X = 0:(ndistinct.Z - 1),
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
  et <- Sys.time()
  time.xz.param <- as.numeric(difftime(et, st, units = "secs"))

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

  ## nonparametric estimated distribution of X|Z (using B-splines)
  st <- Sys.time()
  spline.res.x <- fit.spline(dat = dat, m.knots = m.knots,
                             deg = deg, Boundary.knots = c(0, 1))
  eta1.nonpar <- spline.res.x$dens
  theta.hat.x <- spline.res.x$theta
  knots.x <- spline.res.x$knots
  wts.x <- spline.res.x$wts
  glogf <- spline.res.x$glogf
  et <- Sys.time()
  time.xz.nonpar <- as.numeric(difftime(et, st, units = "secs"))

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
  st <- Sys.time()
  c.params.hat.correct <- t(vapply(
    X = 0:(ndistinct.Z - 1),
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
  et <- Sys.time()
  time.cz.param <- as.numeric(difftime(et, st, units = "secs"))

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

  ## nonparametric estimated distribution of C|Z (using B-splines)
  st <- Sys.time()
  spline.res.c <- fit.spline(dat = mutate(dat, Delta = 1 - Delta),
                             m.knots = m.knots, deg = deg,
                             Boundary.knots = c(0, 1))
  et <- Sys.time()
  time.cz.nonpar <- as.numeric(difftime(et, st, units = "secs"))

  eta2.nonpar <- spline.res.c$dens

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

  ## evaluate CC score at truth and naive root
  #get.Scc(dat = dat,
  #        theta = c(B, log(s2)))
  #get.Scc(dat = dat,
  #        theta = c(naive.lm$coef, log(var(naive.lm$resid))))

  ## complete case
  st <- Sys.time()
  B.cc <- get.root(dat = dat, score = get.Scc,
                   start = c(naive.lm$coef, log(var(naive.lm$resid))))
  et <- Sys.time()
  time.est.cc <- as.numeric(difftime(et, st, units = "secs"))
  st <- Sys.time()
  V.cc <- var.est.sand(dat = datcc, theta = B.cc, args = list(),
                       n = sum(dat$Delta),
                       get.S = get.Scc, return.se = T)
  et <- Sys.time()
  time.var.cc <- as.numeric(difftime(et, st, units = "secs"))

  ## oracle
  st <- Sys.time()
  B.or <- get.root(dat = dat0, score = get.Scc, start = B.cc)
  et <- Sys.time()
  time.est.or <- as.numeric(difftime(et, st, units = "secs"))
  st <- Sys.time()
  V.or <- var.est.sand(dat = dat0, theta = B.or, args = list(),
                       get.S = get.Scc, return.se = T)
  et <- Sys.time()
  time.var.or <- as.numeric(difftime(et, st, units = "secs"))

  ## semiparametric (X|Z, C|Z correct)
  sp.args <- list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                  x.nds = x.nds, c.nds = c.nds,
                  y.nds = y.nds, y.wts = y.wts)

  st <- Sys.time()
  B.sp.11 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.correct,
                                                  "c.wts" = c.wts.correct,
                                                  "eta1" = eta1.correct)))
  et <- Sys.time()
  time.est.sp <- as.numeric(difftime(et, st, units = "secs"))
  st <- Sys.time()
  V.sp.11 <- var.est.sand(dat = dat, theta = B.sp.11,
                          args = append(sp.args, list("x.wts" = x.wts.correct,
                                                      "c.wts" = c.wts.correct,
                                                      "eta1" = eta1.correct)),
                          get.S = get.Seff, return.se = T)
  et <- Sys.time()
  time.var.sp <- as.numeric(difftime(et, st, units = "secs"))

  ## semiparametric (X|Z, C|Z nonparametric)
  B.sp.22 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.nonpar,
                                                  "c.wts" = c.wts.nonpar,
                                                  "eta1" = eta1.nonpar)))
  V.sp.22 <- var.est.sand(dat = dat, theta = B.sp.22,
                          args = append(sp.args, list("x.wts" = x.wts.nonpar,
                                                      "c.wts" = c.wts.nonpar,
                                                      "eta1" = eta1.nonpar)),
                          get.S = get.Seff, return.se = T)

  # return setup parameters, estimates, and standard errors (length 47)
  ret <- c(n = n, q = q, ndistinct.Z = ndistinct.Z, seed = seed,
           B.or = B.or, B.cc = B.cc,
           B.sp.11 = B.sp.11,
           B.sp.22 = B.sp.22,
           V.or = V.or, V.cc = V.cc,
           V.sp.11 = V.sp.11,
           V.sp.22 = V.sp.22,
           time.xz.param = time.xz.param, time.xz.nonpar = time.xz.nonpar,
           time.cz.param = time.cz.param, time.cz.nonpar = time.cz.nonpar,
           time.est.or = time.est.or, time.var.or = time.var.or,
           time.est.cc = time.est.cc, time.var.cc = time.var.cc,
           time.est.sp = time.est.sp, time.var.sp = time.var.sp)

  return(ret)
}


