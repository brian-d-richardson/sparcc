###############################################################################
###############################################################################

# Run Simulation 4

# Anonymous

# 2025-05-02

# Purpose: run simulations with varying number of B-spline knots

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
library(dplyr)
library(numDeriv)
library(zipfR)
load_all()

# simulation parameters ---------------------------------------------------

base.seed <- 10^6 * as.integer(args) # baseline seed (specific to cluster)
len.out <- 16                        # output size

n.sim <- 1                           # number of sims per cluster
n <- 8000                            # sample size
q <- 0.8                             # censoring proportion
m.knots <- c(1, 10)                  # number of distinct Z values

# run simulations ---------------------------------------------------------

## create simulation input
sim.in <- expand.grid(n = n,
                      q = q,
                      m.knots = m.knots,
                      sim.id = 1:n.sim + base.seed)

## run simulations (roughly 10 minutes per replicate)
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    tryCatch(
      expr =
        sim4(
        n = sim.in$n[ii],
        q = sim.in$q[ii],
        m.knots = sim.in$m.knots[ii],
        seed = sim.in$sim.id[ii]),
      error = function(e) c(
        n = sim.in$n[ii],
        q = sim.in$q[ii],
        m.knots = sim.in$m.knots[ii],
        seed = sim.in$sim.id[ii],
        rep(NA, len.out - 3)))

  },
  FUN.VALUE = numeric(len.out)) |>
  t()


## save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim4_data/sd", as.integer(args), ".csv"))

