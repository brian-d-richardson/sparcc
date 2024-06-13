###############################################################################
###############################################################################

# Estimating function developer

# Brian Richardson

# 2024-01-20

# Purpose: develop estimating functions when X and C have gamma distributions

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
library(statmod)
library(devtools)
library(conf)
library(ggplot2)
load_all()

# define parameters -------------------------------------------------------

set.seed(2)
n <- 10000                 # sample size
s2 <- 0.16                 # variance of Y|X,Z
q <- .5                    # censoring proportion
B <- c(1, 2, -1)           # outcome model parameters
x.means <- c(0.5, 1)       # mean of X at levels of Z
x.shape <- 4               # shape parameter for X
c.shape <- 4               # shape parameter for C
mx <- 40                   # number of nodes for X
mc <- 40                   # number of nodes for C
my <- 5                    # number of nodes for Y
specify.x.gamma <- T       # indicator for estimating X as gamma
specify.c.gamma <- T       # indicator for estimating C as gamma

x.rates <- x.shape / x.means  # rate parameters for gamma distribution of X|Z
c.rates <- vapply(
  X = x.rates,                # rate parameters for gamma distribution of C|Z
  FUN.VALUE = 0,
  FUN = function(xr)
    get.c.rate(q = q, x.rate = xr,
               x.shape = x.shape,
               c.shape = c.shape))

# generate data -----------------------------------------------------------

dat.list <- gen.data.gamma(n = n, q = q, B = B, s2 = s2,
                           x.means = x.means,
                           x.shape = x.shape, c.shape = c.shape)
datf <- dat.list$datf          # full data
dat0 <- dat.list$dat0          # oracle data
dat <- dat.list$dat            # observed data
datcc <- dat.list$datcc        # complete case data
zs <- sort(unique(dat$Z))      # unique z values

# plot generated data
ggplot(datf,
       aes(x = X,
           y = Y,
           color = factor(Z))) +
  geom_point() +
  labs(color = "Z")

# estimate nuisance distributions -----------------------------------------

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

# create quadrature rules -------------------------------------------------

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

# evaluate estimating functions -------------------------------------------

# oracle
S0 <- get.Scc(dat = dat0, B = B, ls2 = log(s2),
              args = list(mu = mu, d.mu = d.mu, SF = SF),
              return.sums = F)
assess.ee(S0)

# complete case
Scc <- get.Scc(dat = dat, B = B, ls2 = log(s2),
               args = list(mu = mu, d.mu = d.mu, SF = SF),
               return.sums = F)
assess.ee(Scc)

# MLE
Sml <- get.Sml(dat = dat, B = B, ls2 = log(s2),
             args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                         x.nds = x.nds, x.wts = x.wts),
             return.sums = F)
assess.ee(Sml)

# semiparametric efficient score
Seff <- get.Seff(dat = dat, B = B, ls2 = log(s2),
                 args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy, eta1 = eta1,
                             x.nds = x.nds, x.wts = x.wts,
                             c.nds = c.nds, c.wts = c.wts,
                             y.nds = y.nds, y.wts = y.wts),
                 return.sums = F)
assess.ee(Seff)

# estimate beta -----------------------------------------------------------

# complete case lm to get starting value
naive.lm <- lm(Y ~ W + Z, data = datcc)

# complete case
Bcc <- get.root(dat = dat, score = get.Scc,
                start = c(naive.lm$coef, log(var(naive.lm$resid))))
Bcc

# oracle
B0 <- get.root(dat = dat0, score = get.Scc,
               start = Bcc)
B0

# MLE
Bmle <- get.root(dat = dat, score = get.Sml, start = Bcc,
                 args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                             x.nds = x.nds, x.wts = x.wts))
Bmle

# semiparametric efficient score
Beff <- get.root(dat = dat, score = get.Seff, start = Bcc,
                 args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                             eta1 = eta1,
                             x.nds = x.nds, x.wts = x.wts,
                             c.nds = c.nds, c.wts = c.wts,
                             y.nds = y.nds, y.wts = y.wts))
Beff

# compare estimates
rbind(c(B, log(s2)), B0, Bcc, Bmle, Beff)



