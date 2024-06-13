###############################################################################
###############################################################################

# Nonparametric nuisance estimation developer

# Brian Richardson

# 2024-04-09

# Purpose: estimation procedure with survPresmooth nonparametric estimation

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
library(statmod)
library(devtools)
library(fitdistrplus)
library(ggplot2)
library(survPresmooth)
load_all()

# define parameters -------------------------------------------------------

set.seed(1)
n <- 8000              # sample size
q <- 0.8              # censoring proportion
B <- c(1, 10, 2)      # beta
s2 <- 1               # variance of Y|X,Z
ls2 <- log(s2)
x.thetas <- c(0.5, -0.5)   # theta parameters for X|Z
x.gamma <- 1          # gamma parameter for X|Z
c.gamma <- 2          # gamma parameter for C|Z
mx <- 40                   # number of nodes for X
mc <- 40                   # number of nodes for C
my <- 5                    # number of nodes for Y

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

# bandwidths obtained once using bw.selec = "plug-in" (very slow)
#x.bw <- list(c(0.000842, 0.093202), c(0.0279, 0.2089))
#c.bw <- list(c(0.0668, 0.0408), c(0.000238, 0.106667))

# bandwidths set to reasonable looking values
x.bw <- list(c(0.2, 0.2), c(0.2, 0.2))
c.bw <- list(c(0.2, 0.2), c(0.2, 0.2))

# nonparametrically estimate distribution of X|Z

x.dens.hat <- lapply(
  as.list(zs),
  function(z) survPresmooth::presmooth(
    times = W,
    status = Delta,
    dataset = dat[dat$Z == z,],
    x.est = seq(0, 1, length = 10^3),
    estimand = "f",
    fixed.bw = x.bw[[z + 1]],
    #bw.selec = "plug-in",
    ))

eta1 <- function(x, z) {
  eta <- numeric(length(x))
  for (zz in unique(z)) {
    eta[z == zz] <- approx(
      x = x.dens.hat[[zz + 1]]$x.est,
      y = x.dens.hat[[zz + 1]]$estimate,
      xout = x[z == zz],
      rule = 2
    )$y
  }
  return(eta)
}

# nonparametrically estimate distribution of C|Z

c.dens.hat <- lapply(
  as.list(zs),
  function(z) survPresmooth::presmooth(
    times = dat$W[dat$Z == z],
    status = 1 - dat$Delta[dat$Z == z],
    estimand = "f",
    x.est = seq(0, 1, length = 10^3),
    fixed.bw = c.bw[[z + 1]],
    #bw.selec = "plug-in",
    ))

eta2 <- function(c, z) {
  eta <- numeric(length(c))
  for (zz in unique(z)) {
    eta[z == zz] <- approx(
      x = c.dens.hat[[zz + 1]]$x.est,
      y = c.dens.hat[[zz + 1]]$estimate,
      xout = c[z == zz],
      rule = 2
    )$y
  }
  return(eta)
}

# plot densities ----------------------------------------------------------

# grid to plot X density
x.grid <- seq(0, 1, length = 100)

# plot data
datf %>%
  mutate(e1 = eta1(x = X, z = Z),
         e2 = eta2(c = C, z = Z)) %>%
  ggplot() +
  geom_histogram(aes(x = X,
                     y = after_stat(density)),
                 fill = "blue", alpha = 0.5, bins = n/50) +
  geom_histogram(aes(x = C,
                     y = after_stat(density)),
                 fill = "red", alpha = 0.5, bins = n/50) +
  geom_line(aes(x = X,
                y = e1),
            color = "blue", linewidth = 1) +
  geom_line(aes(x = C,
                y = e2),
            color = "red", linewidth = 1) +
  facet_wrap(~ Z,
             scales = "free", labeller = label_both) +
  labs(x = "X (blue) or C (red)",
       y = "Density") +
  ggtitle("Estimated vs Observed Distributions of X|Z and C|Z")

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
rbind(c(B, log(s2)), B0, Bcc, Bmle, Beff)




