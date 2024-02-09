###############################################################################
###############################################################################

# Simulation 1

# Brian Richardson

# 2024-01-20

# Purpose: run simulations comparing methods under misspecification

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
library(devtools)
#install_github("brian-d-richardson/sparcc"); library(sparcc)
setwd(dirname(getwd()))
load_all()

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
args <- commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 1

# output size
len.out <- 20

# varied parameters
n <- 10000
q <- c(0.4, 0.8)
x.shape <- c(1, 4)
c.shape <- c(1, 4)

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      q = q,
                      x.shape = x.shape,
                      c.shape = c.shape,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    sim1(n = sim.in$n[ii],
         q = sim.in$q[ii],
         x.shape = sim.in$x.shape[ii],
         c.shape = sim.in$c.shape[ii],
         mx = 100,
         mc = 15,
         my = 2,
         seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(len.out)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim1/sd", as.integer(args), ".csv"))




