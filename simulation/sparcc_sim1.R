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
library(statmod)
#install_github("brian-d-richardson/sparcc"); library(sparcc)
#setwd(dirname(getwd()))
load_all()

# simulation parameters ---------------------------------------------------

# baseline seed (specific to cluster)
args <- 1#commandArgs(TRUE)
base.seed <- 10^6 * as.integer(args)

# number of sims per cluster
n.sim <- 1

# output size
len.out <- 20

# fixed parameters
x.shape <- 1.2
c.shape <- 1.2
x.mean <- 0.5
s2 <- 1.1

# varied parameters
n <- 10000
q <- c(0.4, 0.8)
specify.x.gamma <- c(T, F)
specify.c.gamma <- c(T, F)

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      q = q,
                      specify.x.gamma = specify.x.gamma,
                      specify.c.gamma = specify.c.gamma,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    sim1(n = sim.in$n[ii],
         q = sim.in$q[ii],
         x.mean = x.mean,
         s2 = s2,
         x.shape = x.shape,
         c.shape = x.shape,
         specify.x.gamma = sim.in$specify.x.gamma[ii],
         specify.c.gamma = sim.in$specify.c.gamma[ii],
         mx = 100,
         mc = 15,
         my = 3,
         seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(len.out)) |>
  t()

# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim1/sd", as.integer(args), ".csv"))




