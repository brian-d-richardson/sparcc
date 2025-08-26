#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' simulation 4: assess the sparcc estimator with poor B-spline performance
#'
#' @inheritParams get.q.beta
#' @inheritParams gen.data.beta
#'
#' @param m.knots a positive integer, the number of B-spline knots
#' @param seed a positive integer, the seed value for random number generation
#'
#' @importFrom statmod gauss.quad
#' @importFrom zipfR Ibeta
#' @import dplyr
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sim4 <- function(n, q, m.knots, seed) {

  ## for troublehsooting
  #rm(list = ls()); library(devtools); load_all(); library(numDeriv)
  #n = 8000; q = 0.8; seed <- 1; m.knots <- 1

  # define parameters -------------------------------------------------------

  # outcome model parameters
  B <- c(1, 10, 2)
  s2 <- 1

  # nuisance model parameters
  x.thetas <- 0.5 * c(-1, 1)
  x.gamma <- 1
  c.gamma <- 2

  # quadrature rule parameters
  mx <- 40
  mc <- 40
  my <- 5

  # spline estimation parameters
  deg <- 3

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

  ## complete case
  st <- Sys.time()
  B.cc <- get.root(dat = dat, score = get.Scc,
                   start = c(naive.lm$coef, log(var(naive.lm$resid))))
  et <- Sys.time()
  time.est.cc <- as.numeric(difftime(et, st, units = "secs"))

  ## semiparametric (X|Z, C|Z nonparametric)
  sp.args <- list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                  x.nds = x.nds, c.nds = c.nds,
                  y.nds = y.nds, y.wts = y.wts,
                  x.wts = x.wts.nonpar,
                  c.wts = c.wts.nonpar,
                  eta1 = eta1.nonpar)
  st <- Sys.time()
  B.sp.22 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = sp.args)
  et <- Sys.time()
  time.est.sp <- as.numeric(difftime(et, st, units = "secs"))
  st <- Sys.time()
  V.sp.22 <- var.est.sand(dat = dat, theta = B.sp.22,
                          args = sp.args,
                          get.S = get.Seff, return.se = T)
  et <- Sys.time()
  time.var.sp <- as.numeric(difftime(et, st, units = "secs"))

  # return setup parameters, estimates, and standard errors (length 16)
  ret <- c(n = n, q = q, m.knots = m.knots, seed = seed,
           B.sp.22 = B.sp.22,
           V.sp.22 = V.sp.22,
           time.xz.nonpar = time.xz.nonpar,
           time.cz.nonpar = time.cz.nonpar,
           time.est.sp = time.est.sp,
           time.var.sp = time.var.sp)

  return(ret)
}


