###############################################################################
###############################################################################

# Develop Data Generation

# Brian Richardson

# 2024-01-22

# Purpose: develop methods to generate data

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
library(ggplot2)
library(devtools)
library(conf)
load_all()

# define parameters -------------------------------------------------------

n <- 10000            # sample size
q <- 0.8              # censoring proportion
B <- c(1, 2)          # beta
s2 <- 1.1             # variance of Y|X,Z
x.mean <- 0.5         # mean of X

# generate data -----------------------------------------------------------

# both correct
assess.dat(n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
           x.shape = 2, c.shape = 2,
           specify.x.gamma = T, specify.c.gamma = T)

# X incorrect, C correct
assess.dat(n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
           x.shape = 2, c.shape = 2,
           specify.x.gamma = F, specify.c.gamma = T)

# X correct, C incorrect
assess.dat(n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
           x.shape = 2, c.shape = 2,
           specify.x.gamma = T, specify.c.gamma = F)

# X both Incorrect
assess.dat(n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
           x.shape = 2, c.shape = 2,
           specify.x.gamma = F, specify.c.gamma = F)

# assess estimation of etas -----------------------------------------------

n.rep <- 500

assess.eta.est(n.rep = n.rep, n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
               x.shape = 2, c.shape = 2,
               specify.x.gamma = T, specify.c.gamma = T)

