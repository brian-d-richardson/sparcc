###############################################################################
###############################################################################

# Simulation 2

# Brian Richardson

# 2024-01-20

# Purpose: run simulations comparing methods under misspecification

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
library(devtools)
library(statmod)
library(fitdistrplus)
library(dplyr)
#install_github("brian-d-richardson/sparcc"); library(sparcc)
#setwd(dirname(getwd()))
load_all()

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
args <- 1#commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 50

# output size
len.out <- 45

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
x.correct <- T#c(T, F)       # indicator for estimating X as gamma
c.correct <- T#c(T, F)       # indicator for estimating C as gamma

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      q = q,
                      B2 = B2,
                      s2 = s2,
                      x.gamma = x.gamma,
                      c.gamma = c.gamma,
                      x.theta.scale = x.theta.scale,
                      x.correct = x.correct,
                      c.correct = c.correct,
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
        x.correct = sim.in$x.correct[ii],
        c.correct = sim.in$c.correct[ii],
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

