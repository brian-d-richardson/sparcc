###############################################################################
###############################################################################

# Estimating function developer

# Brian Richardson

# 2024-01-20

# Purpose: develop estimating functions

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
library(statmod)
library(devtools)
load_all()

# define parameters -------------------------------------------------------

B <- c(1, 2)               # outcome model parameters
s2 <- 1.1                  # Var(Y|X,Z)
q <- 0.8                   # censoring proportion
n <- 10000                 # sample size
x.mean <- 0.25
x.shape <- 10
c.shape <- 10
x.rate <- x.shape / x.mean # rate parameter for gamma distribution of X
c.rate <- get.c.rate(      # rate parameter for gamma distribution of C
  q = q,
  x.rate = x.rate,
  x.shape = x.shape,
  c.shape = c.shape)
mx <- 100                  # nodes in quadrature grid for X
mc <- 15                   # nodes in quadrature grid for C
my <- 3                    # nodes in quadrature grid for Y

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

# generate data -----------------------------------------------------------

dat.list <- gen.data(n = n, q = q, B = B, s2 = s2,
                     x.rate = x.rate, x.shape = x.shape,
                     c.rate = c.rate, c.shape = c.shape)
datf <- dat.list$datf          # full data
dat0 <- dat.list$dat0          # oracle data
dat <- dat.list$dat            # observed data
datcc <- dat.list$datcc        # complete case data

# estimate nuisance distributions -----------------------------------------

x.rate.hat <- mean(dat$Delta) / mean(dat$W)      # mle for exponential X rate
c.rate.hat <- mean(1 - dat$Delta) / mean(dat$W)  # mle for exponential C rate

# X density
eta1 <- function(x) dexp(x, rate = x.rate.hat)

# C density
eta2 <- function(c) dexp(c, rate = c.rate.hat)

# create quadrature rules -------------------------------------------------

# X quadrature
x.nds <- seq(10^-6, max(datf$X), length = mx)
x.wts <- eta1(x.nds) / sum(eta1(x.nds))

# C quadrature
c.nds <- seq(10^-6, max(datf$C), length = mc)
c.wts <- eta2(c.nds) / sum(eta2(c.nds))

# Y quadrature
gq <- gauss.quad(n = my, kind = "hermite")
y.nds <- gq$nodes
y.wts <- gq$weights

# evaluate estimating functions -------------------------------------------

# oracle
S0 <- get.Scc(dat = dat0, B = B, s2 = s2,
              args = list(mu = mu, d.mu = d.mu, SF = SF),
              return.sums = F)
assess.ee(S0)

# complete case
Scc <- get.Scc(dat = dat, B = B, s2 = s2,
               args = list(mu = mu, d.mu = d.mu, SF = SF),
               return.sums = F)
assess.ee(Scc)

# MLE
Sml <- get.Sml(dat = dat, B = B, s2 = s2,
             args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                         x.nds = x.nds, x.wts = x.wts),
             return.sums = F)
assess.ee(Sml)

# semiparametric efficient score
Seff <- get.Seff(dat = dat, B = B, s2 = s2,
                 args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy, eta1 = eta1,
                             x.nds = x.nds, x.wts = x.wts,
                             c.nds = c.nds, c.wts = c.wts,
                             y.nds = y.nds, y.wts = y.wts),
                 return.sums = F)
assess.ee(Seff)

# estimate beta -----------------------------------------------------------

# oracle
B0 <- get.root(dat = dat0, score = get.Scc, start = c(0, 0, 0))
B0

# complete case
Bcc <- get.root(dat = dat, score = get.Scc, start = c(0, 0, 0))
Bcc

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






