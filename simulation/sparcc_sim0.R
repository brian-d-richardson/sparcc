###############################################################################
###############################################################################

# Simulation 0

# Brian Richardson

# 2024-03-26

# Purpose: simulations comparing methods under misspecification with known etas

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
library(devtools)
library(statmod)
library(conf)
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
len.out <- 25

# parameters
n <- 10000
q <- c(0.6, 0.8)
B2 <- 4
s2 <- 1.1
x.shape <- 1.2
c.shape <- 1.2
x.mean <- 0.25
specify.x.gamma <- T#c(T, F)
specify.c.gamma <- T#c(T, F)

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      q = q,
                      B2 = B2,
                      s2 = s2,
                      x.shape = x.shape,
                      c.shape = c.shape,
                      x.mean = x.mean,
                      specify.x.gamma = specify.x.gamma,
                      specify.c.gamma = specify.c.gamma,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    sim0(n = sim.in$n[ii],
         q = sim.in$q[ii],
         B2 = sim.in$B2[ii],
         s2 = sim.in$s2[ii],
         x.mean = sim.in$x.mean[ii],
         x.shape = sim.in$x.shape[ii],
         c.shape = sim.in$x.shape[ii],
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
          paste0("simulation/sim_data/sim0/sd", as.integer(args), ".csv"))




