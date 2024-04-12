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
len.out <- 29

# parameters
n <- 4000
q <- 0.5
B2 <- 0.5
s2 <- 0.81
x.shape <- 1.25
c.shape <- 1.25
x.mean.scale <- 1
specify.x.gamma <- c(T, F)
specify.c.gamma <- c(T, F)

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      q = q,
                      B2 = B2,
                      s2 = s2,
                      x.shape = x.shape,
                      c.shape = c.shape,
                      x.mean.scale = x.mean.scale,
                      specify.x.gamma = specify.x.gamma,
                      specify.c.gamma = specify.c.gamma,
                      sim.id = 1:n.sim + base.seed)

# run simulations
sim.out <- pbapply::pbvapply(
  X = 1:nrow(sim.in),
  FUN = function(ii) {

    sim1(n = sim.in$n[ii],
         q = sim.in$q[ii],
         B2 = sim.in$B2[ii],
         s2 = sim.in$s2[ii],
         x.mean.scale = sim.in$x.mean.scale[ii],
         x.shape = sim.in$x.shape[ii],
         c.shape = sim.in$c.shape[ii],
         specify.x.gamma = sim.in$specify.x.gamma[ii],
         specify.c.gamma = sim.in$specify.c.gamma[ii],
         mx = 50,
         mc = 15,
         my = 3,
         seed = sim.in$sim.id[ii])

  },
  FUN.VALUE = numeric(len.out)) |>
  t()


# save sim results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim1_data/sd", as.integer(args), ".csv"))




