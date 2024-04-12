###############################################################################
###############################################################################

# Simulation 4

# Brian Richardson

# 2024-04-09

# Purpose: run simulations with nonparametric nuisance estimation

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
library(devtools)
library(statmod)
library(fitdistrplus)
library(dplyr)
library(DALSM)
#setwd(dirname(getwd()))
load_all()

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
args <- 1#commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 1

# output size
len.out <- 43

# parameters
n <- 8000
q <- 0.8
B2 <- 10
s2 <- 1
x.theta.scale <- 0.5
x.gamma <- 1          # gamma parameter for X|Z
c.gamma <- 2          # gamma parameter for C|Z
mx <- 40                   # number of nodes for X
mc <- 40                   # number of nodes for C
my <- 5                    # number of nodes for Y

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      q = q,
                      B2 = B2,
                      s2 = s2,
                      x.gamma = x.gamma,
                      c.gamma = c.gamma,
                      x.theta.scale = x.theta.scale,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    tryCatch(
      expr = sim4(
        n = sim.in$n[ii],
        q = sim.in$q[ii],
        B2 = sim.in$B2[ii],
        s2 = sim.in$s2[ii],
        x.theta.scale = sim.in$x.theta.scale[ii],
        x.gamma = sim.in$x.gamma[ii],
        c.gamma = sim.in$c.gamma[ii],
        mx = mx,
        mc = mc,
        my = my,
        seed = sim.in$sim.id[ii]),
      error = function(e) c(
        n = sim.in$n[ii],
        q = sim.in$q[ii],
        B2 = sim.in$B2[ii],
        s2 = sim.in$s2[ii],
        rep(NA, len.out - 4)))

  },
  FUN.VALUE = numeric(len.out)) |>
  t()


# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim4_data/sd", as.integer(args), ".csv"))

