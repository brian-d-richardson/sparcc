###############################################################################
###############################################################################

# Run Simulation 3

# Anonymous

# 2025-04-30

# Purpose: run simulations with high-dimensional Z

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
source("simulation/sim3_function.R")

# simulation parameters ---------------------------------------------------

base.seed <- 10^6 * as.integer(args) # baseline seed (specific to cluster)
len.out <- 46                        # output size

n.sim <- 2                          # number of sims per cluster
n <- 8000                            # sample size
q <- 0.8                             # censoring proportion
ndistinct.Z <- 8                     # number of distinct Z values

# run simulations ---------------------------------------------------------

## create simulation input
sim.in <- expand.grid(n = n,
                      q = q,
                      ndistinct.Z = ndistinct.Z,
                      sim.id = 1:n.sim + base.seed)

## run simulations (roughly 10 minutes per replicate)
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    tryCatch(
      expr =
        sim3(
        n = sim.in$n[ii],
        q = sim.in$q[ii],
        ndistinct.Z = sim.in$ndistinct.Z[ii],
        seed = sim.in$sim.id[ii]),
      error = function(e) c(
        n = sim.in$n[ii],
        q = sim.in$q[ii],
        ndistinct.Z = sim.in$ndistinct.Z[ii],
        seed = sim.in$sim.id[ii],
        rep(NA, len.out - 4)))

  },
  FUN.VALUE = numeric(len.out)) |>
  t()


## save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim3_data/sd", as.integer(args), ".csv"))

