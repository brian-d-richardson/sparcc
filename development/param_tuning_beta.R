###############################################################################
###############################################################################

# Estimating function developer

# Brian Richardson

# 2024-04-02

# Purpose: tune parameters for efficient score function when X, C are beta

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
library(statmod)
library(devtools)
library(fitdistrplus)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pbapply)
load_all()

# define parameters -------------------------------------------------------

set.seed(2)
n <- 8000             # sample size
q <- 0.5              # censoring proportion
B <- c(1, 10, 2)      # beta
s2 <- 4               # variance of Y|X,Z
x.thetas <- c(0.5, -0.5)   # theta parameters for X|Z
x.gamma <- 1          # gamma parameter for X|Z
c.gamma <- 4          # gamma parameter for C|Z
x.correct <- T       # indicator for estimating X as gamma
c.correct <- T       # indicator for estimating C as gamma

# generate data -----------------------------------------------------------

dat.list <- gen.data.beta(
  n = n, q = q, B = B, s2 = s2,
  x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma)

datf <- dat.list$datf          # full data
dat0 <- dat.list$dat0          # oracle data
dat <- dat.list$dat            # observed data
datcc <- dat.list$datcc        # complete case data
zs <- sort(unique(dat$Z))      # unique z values

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

# search over grid of mx, mc, my ------------------------------------------

run.search <- T

# grid of possible values for mx and mc
mx.grid <- seq(15, 55, by = 5)
mc.grid <- seq(15, 55, by = 5)
my.grid <- seq(2, 6, by = 2)
search.in <- expand.grid(mx = mx.grid,
                         mc = mc.grid,
                         my = my.grid)

# store Seff and computation time for each mx, mc, my (takes ~ 4 min to run)
if (run.search) {

  search.out <- pbvapply(
    X = 1:nrow(search.in),
    FUN.VALUE = numeric(8),
    FUN = function(ii) {

      # X|Z quadrature
      x.nds <- vapply(
        X = zs,
        FUN.VALUE = numeric(search.in$mx[ii]),
        FUN = function(z) seq(1E-6, 1-1E-6, length = search.in$mx[ii]))

      x.wts <- vapply(
        X = 1:length(zs),
        FUN.VALUE = numeric(search.in$mx[ii]),
        FUN = function(i) eta1(x.nds[,i], zs[i]) / sum(eta1(x.nds[,i], zs[i])))

      # C|Z quadrature
      c.nds <- vapply(
        X = zs,
        FUN.VALUE = numeric(search.in$mc[ii]),
        FUN = function(z) seq(1E-6, 1-1E-6, length = search.in$mc[ii]))

      c.wts <- vapply(
        X = 1:length(zs),
        FUN.VALUE = numeric(search.in$mc[ii]),
        FUN = function(i) eta2(c.nds[,i], zs[i]) / sum(eta2(c.nds[,i], zs[i])))

      # Y quadrature
      gq <- gauss.quad(n = search.in$my[ii], kind = "hermite")
      y.nds <- gq$nodes
      y.wts <- gq$weights

      # evaluate efficient score
      st <- Sys.time()
      Seff <- get.Seff(dat = dat, B = B, ls2 = log(s2),
                       args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy, eta1 = eta1,
                                   x.nds = x.nds, x.wts = x.wts,
                                   c.nds = c.nds, c.wts = c.wts,
                                   y.nds = y.nds, y.wts = y.wts),
                       return.sums = T)
      et <- Sys.time()

      return(c(mx = search.in$mx[ii],
               mc = search.in$mc[ii],
               my = search.in$my[ii],
               Seff = Seff,
               Time = et - st))
    }) %>%
    t() %>%
    as.data.frame()

  write.csv(search.out, "development/dev_data/m_tuning_beta.csv", row.names = F)
}

# plot results ------------------------------------------------------------

# load results
search.out <- read.csv("development/dev_data/m_tuning_beta.csv")
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


