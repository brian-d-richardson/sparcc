###############################################################################
###############################################################################

# Run Simulation 2

# Anonymous

# 2024-04-04

# Purpose: run simulations comparing methods over sequence of q

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())

# indicator for running on cluster
on.cluster <- F
if (on.cluster) {
  setwd(dirname(getwd()))
  args <- commandArgs(TRUE)          # cluster id
} else {
  args <- 1
}

library(devtools)
library(statmod)
library(fitdistrplus)
library(dplyr)
load_all()

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 100

# output size
len.out <- 27

# parameters
n <- c(200, 800, 8000)
q <- c(0.1, 0.2)#c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
       #0.65, 0.7, 0.75, 0.8, 0.85, 0.9)
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
      expr = sim2(
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
          paste0("simulation/sim2_data/sd", as.integer(args), ".csv"))

