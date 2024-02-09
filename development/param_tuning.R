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
setwd(dirname(getwd()))
library(statmod)
library(ggplot2)
library(tidyverse)
library(pbapply)
source("R/int_eq.R"); source("R/est_funs.R"); source("R/misc_helpers.R")

# define parameters -------------------------------------------------------

n <- 10000                 # sample size
q <- .8                    # censoring proportion
B <- c(1, 2)               # outcome model parameters
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
fy <- function(y, x, B) dnorm(x = y, mean = mu(x, B), sd = 1)

# full data score vector
SBF <- function(y, x, B) {
  (y - mu(x, B)) *
   d.mu(x, B)
}

# generate data -----------------------------------------------------------

X <- rexp(n, x.rate)                         # censored covariate
C <- rexp(n, c.rate)                         # censoring time
W <- ifelse(X <= C, X, C)                    # observed covariate
Delta <- ifelse(X <= C, 1, 0)                # uncensored indicator
Y <- rnorm(n, cbind(1, X) %*% B, sd = 1)     # outcome
dat0 <- data.frame(Y, W = X, Delta = 1)      # oracle data
dat <- data.frame(Y, W, Delta)               # observed data
datcc <- dat[Delta == 1,]                    # complete case data


# search over grid of mx, mc, my ------------------------------------------

# grid of possible values for mx and mc
mx.grid <- seq(20, 80, by = 10)
mc.grid <- seq(10, 20, by = 5)
my.grid <- seq(2, 5, by = 1)
search.in <- expand.grid(mx = mx.grid,
                         mc = mc.grid,
                         my = my.grid)

# store Seff and computation time for each mx, mc, my (takes ~ 2 min to run)
search.out <- pbvapply(
  X = 1:nrow(search.in),
  FUN.VALUE = numeric(6),
  FUN = function(ii) {

    # X quadrature
    x.nds <- seq(10^-6, max(X), length = search.in$mx[ii])
    x.wts <- eta1(x.nds) / sum(eta1(x.nds))

    # C quadrature
    c.nds <- seq(10^-6, max(C), length = search.in$mc[ii])
    c.wts <- eta2(c.nds) / sum(eta2(c.nds))

    # Y quadrature
    gq <- gauss.quad(n = search.in$my[ii], kind = "hermite")
    y.nds <- gq$nodes
    y.wts <- gq$weights

    # evaluate efficient score
    st <- Sys.time()
    Seff <- get.Seff(
      dat = dat, B = B,
      args = list(mu = mu, d.mu = d.mu, SBF = SBF, fy = fy, eta1 = eta1,
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

#write.csv(search.out, "sim_data/param_tuning/mc_mc_res.csv", row.names = F)

# plot results ------------------------------------------------------------

# load results
#search.out <- read.csv("sim_data/param_tuning/mc_mc_res.csv")
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


