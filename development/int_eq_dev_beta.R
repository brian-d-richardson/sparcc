###############################################################################
###############################################################################

# Develop integral equation solver

# Brian Richardson

# 2024-04-02

# Purpose: develop integral equation solver

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
mx <- 40                   # number of nodes for X
mc <- 40                   # number of nodes for C
my <- 5                    # number of nodes for Y
x.correct <- F       # indicator for estimating X as gamma
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
SF <- function(y, x, z, B, s2) {
  cbind((y - mu(x, z, B)) * d.mu(x, z, B),
        (y - mu(x, z, B)) ^ 2 - s2)
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
  ggtitle("Oracle Data")

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

# margin for X and C nodes
margin <- 1E-6

# X|Z quadrature
x.nds <- vapply(
  X = zs,
  FUN.VALUE = numeric(mx),
  FUN = function(z) seq(margin, 1 - margin, length = mx))

x.wts <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mx),
  FUN = function(i) eta1(x.nds[,i], zs[i]) / sum(eta1(x.nds[,i], zs[i])))

# C|Z quadrature
c.nds <- vapply(
  X = zs,
  FUN.VALUE = numeric(mc),
  FUN = function(z) seq(margin, 1 - margin, length = mc))

c.wts <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mc),
  FUN = function(i) eta2(c.nds[,i], zs[i]) / sum(eta2(c.nds[,i], zs[i])))

# Y quadrature
gq <- gauss.quad(n = my, kind = "hermite")
y.nds <- gq$nodes
y.wts <- gq$weights

# solve integral equation -------------------------------------------------

# solve integral equation for a at node points, for each unique Z value
a <- get.a(
  B = B, s2 = s2, mu = mu, d.mu = d.mu,
  fy = fy, SF = SF, zs = zs,
  x.nds = x.nds, x.wts = x.wts,
  c.nds = c.nds, c.wts = c.wts,
  y.nds = y.nds, y.wts = y.wts)

# check that a has mean zero
colMeans(a[[1]] * x.wts[,1])
colMeans(a[[2]] * x.wts[,2])

# plot a
adat <- data.frame(rbind(
  cbind(zs[1], x.nds[,1], a[[1]], a[[1]]*x.wts[,1]),
  cbind(zs[2], x.nds[,2], a[[2]], a[[2]]*x.wts[,2]))) |>
  `colnames<-`(c("z", "x",
                 "a1", "a2", "a3", "a4",
                 "a.eta1", "a.eta2", "a.eta3", "a.eta4")) |>
  pivot_longer(cols = !c(x, z)) |>
  mutate(comp = gsub("\\D", x = name, replacement = ""),
         func = gsub("\\d", x = name, replacement = ""))

ggplot(adat,
       aes(x = x,
           y = value,
           color = comp)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  facet_wrap(z ~ func,
             scales = "free",
             labeller = label_both)
