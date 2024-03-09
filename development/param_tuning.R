###############################################################################
###############################################################################

# Estimating function developer

# Brian Richardson

# 2024-01-24

# Purpose: tune parameters for efficient score function

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
library(statmod)
library(ggplot2)
library(tidyverse)
library(pbapply)
library(devtools)
load_all()

# define parameters -------------------------------------------------------

n <- 8000                  # sample size
q <- 0.8                   # censoring proportion
B <- c(1, .2)              # outcome model parameters
s2 <- 1.1                  # Var(Y|X,Z)
x.mean <- 0.5
x.shape <- 1.2
c.shape <- 2
x.rate <- x.shape / x.mean # rate parameter for gamma distribution of X
c.rate <- get.c.rate(      # rate parameter for gamma distribution of C
  q = q,
  x.rate = x.rate,
  x.shape = x.shape,
  c.shape = c.shape)
specify.x.gamma <- T
specify.c.gamma <- T

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

## generate data
dat.list <- gen.data(n = n, q = q, B = B, s2 = s2,
                     x.rate = x.rate, x.shape = x.shape,
                     c.rate = c.rate, c.shape = c.shape)
datf <- dat.list$datf          # full data
dat0 <- dat.list$dat0          # oracle data
dat <- dat.list$dat            # observed data
datcc <- dat.list$datcc        # complete case data

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

x.upper <- max(datf$X)  # oracle max
c.upper <- max(datf$C)

# search over grid of mx, mc, my ------------------------------------------

# grid of possible values for mx and mc
mx.grid <- seq(40, 120, by = 10)
mc.grid <- seq(10, 20, by = 5)
my.grid <- seq(2, 4, by = 1)
search.in <- expand.grid(mx = mx.grid,
                         mc = mc.grid,
                         my = my.grid)

# store Seff and computation time for each mx, mc, my (takes ~ 4 min to run)
search.out <- pbvapply(
  X = 1:nrow(search.in),
  FUN.VALUE = numeric(7),
  FUN = function(ii) {

    # X quadrature
    x.nds <- seq(10^-6, max(datf$X), length = search.in$mx[ii])
    x.wts <- eta1(x.nds) / sum(eta1(x.nds))

    # C quadrature
    c.nds <- seq(10^-6, max(datf$C), length = search.in$mc[ii])
    c.wts <- eta2(c.nds) / sum(eta2(c.nds))

    # Y quadrature
    gq <- gauss.quad(n = search.in$my[ii], kind = "hermite")
    y.nds <- gq$nodes
    y.wts <- gq$weights

    # evaluate efficient score
    st <- Sys.time()
    Seff <- get.Seff(
      dat = dat, B = B, s2 = s2,
      args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy, eta1 = eta1,
                  x.nds = x.nds, x.wts = x.wts,
                  c.nds = c.nds, c.wts = c.wts,
                  y.nds = y.nds, y.wts = y.wts))
    et <- Sys.time()

    return(c(mx = search.in$mx[ii],
             mc = search.in$mc[ii],
             my = search.in$my[ii],
             Seff = Seff,
             Time = et - st))
  }) %>%
  t() %>%
  as.data.frame()

#write.csv(search.out, "development/dev_data/m_tuning.csv", row.names = F)

# plot results ------------------------------------------------------------

# load results
#search.out <- read.csv("development/dev_data/m_tuning.csv")
search.out.long <- search.out %>%
  pivot_longer(cols = c(Seff1, Seff2, Time)) %>%
  mutate(mc = factor(mc),
         mx = factor(mx),
         my = factor(my))

mc.labs <- paste0("mc = ", mc.grid); names(mc.labs) <- mc.grid
my.labs <- paste0("my = ", my.grid); names(my.labs) <- my.grid

# facet by Y
ggplot(data = search.out.long,
       aes(x = mx,
           y = value,
           color = mc,
           group = mc)) +
  geom_line() +
  facet_grid(name ~ my,
             scales = "free",
             labeller = labeller(my = my.labs)) +
  labs(y = "") +
  ggtitle("Score Values and Computation Time by Node Counts mx and mc")

# facet by C
ggplot(data = search.out.long,
       aes(x = mx,
           y = value,
           color = my,
           group = my)) +
  geom_line() +
  facet_grid(name ~ mc,
             scales = "free",
             labeller = labeller(mc = mc.labs)) +
  labs(y = "") +
  ggtitle("Score Values and Computation Time by Node Counts mx and mc")


