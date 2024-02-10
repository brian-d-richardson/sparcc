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
install_github("brian-d-richardson/sparcc")
library(sparcc)
library(statmod)

# define parameters -------------------------------------------------------

n <- 10000                 # sample size
s2 <- 1.1                  # variance of Y|X,Z
q <- .6                    # censoring proportion
B <- c(1, 2)               # outcome model parameters
mx <- 21                   # number of nodes for X
mc <- 21                   # number of nodes for C
my <- 11                   # number of nodes for Y
x.rate <- 4                # rate parameter for gamma distribution of X
c.rate <- x.rate*q/(1-q)   # rate parameter for exp distribution of C

# X density
eta1 <- function(x) dexp(x, rate = x.rate)

# C density
eta2 <- function(c) dexp(c, rate = c.rate)

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

X <- rexp(n, x.rate)                            # censored covariate
C <- rexp(n, c.rate)                            # censoring time
W <- ifelse(X <= C, X, C)                       # observed covariate
Delta <- ifelse(X <= C, 1, 0)                   # uncensored indicator
Y <- rnorm(n, cbind(1, X) %*% B, sd = sqrt(s2)) # outcome
dat0 <- data.frame(Y, W = X, Delta = 1)         # oracle data
dat <- data.frame(Y, W, Delta)                  # observed data
datcc <- dat[Delta == 1,]                       # complete case data

# create quadrature rules -------------------------------------------------

# X quadrature
x.nds <- seq(10^-6, max(X), length = mx)
x.wts <- eta1(x.nds) / sum(eta1(x.nds))

# C quadrature
c.nds <- seq(10^-6, max(C), length = mc)
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

# ROOT NOT FOUND WHEN S2 UNKOWN??

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
                             x.nds = x.nds, x.wts = x.wts,
                             c.nds = c.nds, c.wts = c.wts,
                             y.nds = y.nds, y.wts = y.wts))
Beff

# compare estimates
rbind(c(B, log(s2)), B0, Bcc, Bmle, Beff)






