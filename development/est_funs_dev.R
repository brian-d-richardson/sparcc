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
library(conf)
load_all()

# define parameters -------------------------------------------------------

set.seed(1)
B <- c(1, 0.5)            # outcome model parameters
s2 <- 0.25                # Var(Y|X,Z)
q <- 0.8                  # censoring proportion
n <- 10000                # sample size
x.mean <- 1
x.shape <- 1.1
c.shape <- 1.1
x.rate <- x.shape / x.mean # rate parameter for gamma distribution of X
c.rate <- get.c.rate(      # rate parameter for gamma distribution of C
  q = q,
  x.rate = x.rate,
  x.shape = x.shape,
  c.shape = c.shape)
mx <- 100                  # nodes in quadrature grid for X
mc <- 15                   # nodes in quadrature grid for C
my <- 3                    # nodes in quadrature grid for Y
specify.x.gamma <- T       # indicator for estimating X as gamma
specify.c.gamma <- T       # indicator for estimating C as gamma

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

plot(datf$X, datf$Y)
plot(datcc$W, datcc$Y)
max(dat$W[dat$Delta == 1]); max(datf$X)

# estimate nuisance distributions -----------------------------------------

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


# create quadrature rules -------------------------------------------------

x.upper <- max(datf$X)  # oracle max
c.upper <- max(datf$C)

#x.upper <- sum(1/(1:n)) / x.rate.hat # estimated expected max
#c.upper <- sum(1/(1:n)) / c.rate.hat
#x.upper <- qexp((n-1)/n, rate = x.rate.hat)  # estimated (n-1)/n quantile
#c.upper <- qexp((n-1)/n, rate = c.rate.hat)

# X quadrature
x.nds <- seq(10^-5, x.upper, length = mx)
x.wts <- eta1(x.nds) / sum(eta1(x.nds))

# C quadrature
c.nds <- seq(10^-5, c.upper, length = mc)
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

# complete case lm to get starting value
naive.lm <- lm(Y ~ W, data = datcc)

# complete case
Bcc <- get.root(dat = dat, score = get.Scc,
                start = c(naive.lm$coef, log(var(naive.lm$resid))))
Bcc

# oracle
B0 <- get.root.notrycatch(dat = dat0, score = get.Scc,
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






