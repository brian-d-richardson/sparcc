###############################################################################
###############################################################################

# Run Simulation 1

# Brian Richardson

# 2024-04-12

# Purpose: run simulations comparing methods under misspecification

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
library(numDeriv)
library(zipfR)
load_all()

# simulation parameters ---------------------------------------------------

base.seed <- 10^6 * as.integer(args) # baseline seed (specific to cluster)
len.out <- 85                        # output size

n.sim <- 100                         # number of sims per cluster
n <- 8000                            # sample size
q <- c(0.4, 0.8)                     # censoring proportion

# run simulations ---------------------------------------------------------

## create simulation input
sim.in <- expand.grid(n = n,
                      q = q,
                      sim.id = 1:n.sim + base.seed)

## run simulations (roughly 6 minutes per replicate)
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    tryCatch(
      expr =
        sim1(
        n = sim.in$n[ii],
        q = sim.in$q[ii],
        seed = sim.in$sim.id[ii]),
      error = function(e) c(
        n = sim.in$n[ii],
        q = sim.in$q[ii],
        seed = sim.in$sim.id[ii],
        rep(NA, len.out - 3)))

  },
  FUN.VALUE = numeric(len.out)) |>
  t()


## save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim1_data/sd", as.integer(args), ".csv"))

