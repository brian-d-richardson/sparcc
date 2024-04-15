###############################################################################
###############################################################################

# Estimating function developer

# Brian Richardson

# 2024-01-20

# Purpose: develop estimating functions where X|Z, C|Z ~ Beta

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
library(statmod)
library(devtools)
library(fitdistrplus)
library(ggplot2)
load_all()

# define parameters -------------------------------------------------------

set.seed(2)
n <- 8000             # sample size
q <- 0.8              # censoring proportion
B <- c(1, 10, 2)      # beta
s2 <- 2               # variance of Y|X,Z
ls2 <- log(2)
x.thetas <- c(0.5, -0.5)   # theta parameters for X|Z
x.gamma <- 1          # gamma parameter for X|Z
c.gamma <- 2          # gamma parameter for C|Z
mx <- 40                   # number of nodes for X
mc <- 40                   # number of nodes for C
my <- 5                    # number of nodes for Y
x.correct <- T       # indicator for estimating X as gamma
c.correct <- T       # indicator for estimating C as gamma

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
        ((y - mu(x, z, B)) ^ 2 - exp(ls2)) * exp(ls2))
}

# generate data -----------------------------------------------------------

dat.list <- gen.data.beta(
  n = n, q = q, B = B, s2 = s2,
  x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma)

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
  labs(color = "Z") +
  ggtitle("Oracle Data")

ggplot(datcc,
       aes(x = W,
           y = Y,
           color = factor(Z))) +
  geom_point() +
  labs(color = "Z") +
  ggtitle("Complete Case Data")

# estimate nuisance distributions -----------------------------------------

# estimate distribution of X|Z
if (x.correct == T) {

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

} else {

  # misspecify: estimate marginal beta distribution
  est <- dat %>%
    mutate(left = W,
           right = ifelse(Delta == 1, W, NA)) %>%
    dplyr::select(left, right) %>%
    fitdistrplus::fitdistcens(distr = "beta")
  x.params.hat <- est$estimate

  # define estimated X|Z density
  eta1 <- function(x, z) {
    dbeta(x = x,
          shape1 = x.params.hat["shape1"],
          shape2 = x.params.hat["shape2"])
  }
}

# estimate distribution of C|Z
if (c.correct) {

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

} else {

  # misspecify: estimate marginal beta distribution
  est <- dat %>%
    mutate(left = W,
           right = ifelse(Delta == 0, W, NA)) %>%
    dplyr::select(left, right) %>%
    fitdistrplus::fitdistcens(distr = "beta")
  c.params.hat <- est$estimate

  # define estimated X|Z density
  eta2 <- function(c, z) {
    dbeta(x = c,
          shape1 = c.params.hat["shape1"],
          shape2 = c.params.hat["shape2"])
  }
}


# create quadrature rules -------------------------------------------------

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

# evaluate estimating functions -------------------------------------------

# oracle
S0 <- get.Scc(dat = dat0, B = B, ls2 = ls2,
              args = list(mu = mu, d.mu = d.mu, SF = SF),
              return.sums = F)
assess.ee(S0)

# complete case
Scc <- get.Scc(dat = dat, B = B, ls2 = ls2,
               args = list(mu = mu, d.mu = d.mu, SF = SF),
               return.sums = F)
assess.ee(Scc)

# MLE
Sml <- get.Sml(dat = dat, B = B, ls2 = ls2,
             args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                         x.nds = x.nds, x.wts = x.wts),
             return.sums = F)
assess.ee(Sml)

# semiparametric efficient score
Seff <- get.Seff(dat = dat, B = B, ls2 = ls2,
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
rbind(c(B, ls2, B0, Bcc, Bmle, Beff))

round(1E7 * (cbind(B0, Bcc, Bmle, Beff) - c(B, log(s2))) ^ 2)



