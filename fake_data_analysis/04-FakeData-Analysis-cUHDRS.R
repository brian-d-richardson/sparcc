###############################################################################
###############################################################################

# Fake ENROLL-HD Data Analysis with cUHDRS Outcome

# Brian Richardson

# 2025-08-26

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

rm(list = ls())
library(dplyr)
library(ggplot2)
library(devtools)
library(fitdistrplus)
library(MASS)
library(tictoc)
library(statmod)
library(here)
setwd(here())
load_all(dirname(getwd()))

# load data ---------------------------------------------------------------

datf <- read.csv("derived-data/datf_cuhdrs.csv")
dat <- read.csv("derived-data/dat_cuhdrs.csv")

# estimate nuisance distributions -----------------------------------------

zs <- sort(unique(dat$Z)) # unique z values
mx <- 50                  # grid size for X
mc <- 50                  # grid size for C
my <- 5                   # grid size for Y
m.knots <- 5              # number of knots in spline models
Boundary.knots <- c(0, 1) # boundary knots for spline models
deg <- 3                  # degree of spline models

## X|Z quadrature nodes
x.nds <- vapply(
  X = zs,
  FUN.VALUE = numeric(mx),
  FUN = function(z) seq(1E-6, 1-1E-6, length = mx))

## X|Z modeled as conditional gamma
x.params.hat <- t(vapply(
  X = 0:1,
  FUN.VALUE = numeric(2),
  FUN = function(z) {

    est <- dat %>%
      filter(Z == z) %>%
      mutate(left = W,
             right = ifelse(Delta == 1, W, NA)) %>%
      dplyr::select(left, right) %>%
      as.data.frame() %>%
      fitdistcens(distr = "beta")

    return(est$estimate)
  }))

eta1.param <- function(x, z) {
  dbeta(x = x,
        shape1 = x.params.hat[z + 1, "shape1"],
        shape2 = x.params.hat[z + 1, "shape2"])
}

x.wts.param <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mx),
  FUN = function(i) eta1.param(x.nds[,i], zs[i]) /
    sum(eta1.param(x.nds[,i], zs[i])))

## nonparametric estimated distribution of X|Z (using B-splines)
tic("Spline fit X")
spline.res.x <- fit.spline(dat = dat, m.knots = m.knots,
                           Boundary.knots = Boundary.knots, deg = deg)
toc()
eta1.nonpar <- spline.res.x$dens
theta.hat.x <- spline.res.x$theta
knots.x <- spline.res.x$knots
wts.x <- spline.res.x$wts
glogf <- spline.res.x$glogf

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

## C|Z modeled as conditional beta
c.params.hat <- t(vapply(
  X = 0:1,
  FUN.VALUE = numeric(2),
  FUN = function(z) {

    # NOTE: C is available for all subjects; no censoring
    est <- datf %>%
      filter(Z == z) %>%
      dplyr::select(C) %>%
      unlist() %>%
      fitdist(distr = "beta")

    return(est$estimate)
  }))

eta2.param <- function(c, z) {
  dbeta(x = c,
        shape1 = c.params.hat[z + 1, "shape1"],
        shape2 = c.params.hat[z + 1, "shape2"])
}

c.wts.param <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mc),
  FUN = function(i) eta2.param(c.nds[,i], zs[i]) /
    sum(eta2.param(c.nds[,i], zs[i])))

## nonparametric estimated distribution of C|Z (using B-splines)
tic("Spline fit C")
spline.res.c <- fit.spline(dat = mutate(dat, Delta = 1 - Delta),
                           m.knots = m.knots, Boundary.knots = Boundary.knots,
                           deg = deg)
toc()
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

# analysis ----------------------------------------------------------------

## complete case linear model to get starting value
cc.lm <- lm(Y ~ W * Z, data = filter(dat, Delta == 1))
summary(cc.lm)

## complete case
tic("CC estimate")
B.cc <- get.root(dat = dat, score = get.Scc,
                 args = list(xz.interaction = T),
                 start = c(cc.lm$coef, log(var(cc.lm$resid))))
toc()
tic("CC variance")
V.cc <- var.est.sand(dat = dat, theta = B.cc,
                     n = sum(dat$Delta),
                     list(xz.interaction = T),
                     get.S = get.Scc, return.se = F)
toc()
saveRDS(list(B = B.cc, V = V.cc),
        file = "derived-data/cUHDRS-results/cc_res")

cbind(B.cc, sqrt(diag(V.cc)))
summary(cc.lm)$coefficients[,1:2]

## MLE
mle.args <- list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                 x.nds = x.nds, x.wts = x.wts.param,
                 xz.interaction = T)
tic("MLE estimate")
B.ml <- get.root(dat = dat, score = get.Sml, start = B.cc,
                 args = mle.args)
toc()
tic("MLE variance")
V.ml <- var.est.sand(
  dat = dat,
  get.S = function(dat, theta, args, return.sums = F) {

    # extract X|Z parameters
    alpha <- matrix(tail(theta, -5), ncol = length(zs))

    # define estimated X|Z density
    eta1 <- function(x, z) {
      dbeta(x = x,
            shape1 = alpha[z + 1],
            shape2 = alpha[z + 3])
    }

    # create quadrature nodes
    x.wts <- vapply(
      X = 1:length(zs),
      FUN.VALUE = numeric(mx),
      FUN = function(i) eta1(x.nds[,i], zs[i]) / sum(eta1(x.nds[,i], zs[i])))

    args <- list(x.nds = x.nds, x.wts = x.wts, xz.interaction = T)

    # stack estimating equations
    S <- cbind(get.Sml(dat = dat, theta = head(theta, 5),
                       args = args, return.sums = F),
               d.log.fxz(dat = dat, theta = alpha,
                         args = args, return.sums = F))

    if (return.sums) {
      return(colSums(S))
    } else {
      return(S)
    }
  },
  theta = c(B.ml, x.params.hat),
  args = list(),
  return.se = F)
toc()
saveRDS(list(B = B.ml, V = V.ml),
        file = "derived-data/cUHDRS-results/ml_res")

## semiparametric (X|Z, C|Z nonparametric)
sp.nonpar.args <- list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                       x.nds = x.nds, c.nds = c.nds,
                       y.nds = y.nds, y.wts = y.wts,
                       x.wts = x.wts.nonpar,
                       c.wts = c.wts.nonpar,
                       eta1 = eta1.nonpar,
                       xz.interaction = T)
tic("SP (nonpar) estimate")
B.sp.nonpar <- get.root(dat = dat, score = get.Seff,
                        start = B.cc,
                        args = sp.nonpar.args)
toc()
tic("SP (nonpar) variance")
V.sp.nonpar <- var.est.sand(dat = dat, theta = B.sp.nonpar,
                            args = sp.nonpar.args,
                            get.S = get.Seff, return.se = F)
toc()
saveRDS(list(B = B.sp.nonpar, V = V.sp.nonpar),
        file = "derived-data/cUHDRS-results/sp_nonpar_res")

## semiparametric (X|Z, C|Z parametric)
sp.param.args <- list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                      x.nds = x.nds, c.nds = c.nds,
                      y.nds = y.nds, y.wts = y.wts,
                      x.wts = x.wts.param,
                      c.wts = c.wts.param,
                      eta1 = eta1.param,
                      xz.interaction = T)
tic("SP (param) estimate")
B.sp.param <- get.root(dat = dat, score = get.Seff,
                       start = B.sp.nonpar,
                       args = sp.param.args)
toc()
tic("SP (param) variance")
V.sp.param <- var.est.sand(dat = dat, theta = B.sp.param,
                           args = sp.param.args,
                           get.S = get.Seff, return.se = F)
toc()
saveRDS(list(B = B.sp.param, V = V.sp.param),
        file = "derived-data/cUHDRS-results/sp_param_res")

## results
res <- data.frame(
  est = c(B.cc, B.ml, B.sp.param, B.sp.nonpar),
  ste = c(sqrt(diag(V.cc)), sqrt(diag(V.ml)[1:5]),
          sqrt(diag(V.sp.param)), sqrt(diag(V.sp.nonpar))),
  param = rep(1:5, times = 4),
  method = rep(c("CC", "MLE", "Semipar. (Param.)", "Semipar. (Nonpar.)"),
               each = 5)) %>%
  mutate(ci.lower = est - qnorm(0.975) * ste,
         ci.upper = est + qnorm(0.975) * ste)

## save results
write.csv(res, "derived-data/cUHDRS-results/res.csv", row.names = F)






