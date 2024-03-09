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
len.out <- 25

# parameters
n <- 10000
q <- c(0.65, 0.7, 0.75)
B2 <- 0.2
s2 <- 0.09
shape <- c(1.1, 1.15, 1.2)
#x.shape <- 1.1
#c.shape <- 1.1
x.mean <- 0.5
specify.x.gamma <- c(T, F)
specify.c.gamma <- c(T, F)

# run simulations ---------------------------------------------------------

# create simulation input
sim.in <- expand.grid(n = n,
                      q = q,
                      B2 = B2,
                      s2 = s2,
                      shape = shape,
                      #x.shape = x.shape,
                      #c.shape = c.shape,
                      x.mean = x.mean,
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
         x.mean = sim.in$x.mean[ii],
         x.shape = sim.in$shape[ii],
         c.shape = sim.in$shape[ii],
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




