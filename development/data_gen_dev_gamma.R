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
B <- c(1, 2, -1)      # beta
s2 <- 1.1             # variance of Y|X,Z
x.means <- c(0.5, 1)  # means of X at Z = 0, 1

# generate data -----------------------------------------------------------

# both correct
assess.dat.gamma(n = n, q = q, B = B, s2 = s2, x.means = x.means,
           x.shape = 2, c.shape = 2,
           specify.x.gamma = T, specify.c.gamma = T)

# X incorrect, C correct
assess.dat.gamma(n = n, q = q, B = B, s2 = s2, x.means = x.means,
                 x.shape = 2, c.shape = 2,
                 specify.x.gamma = F, specify.c.gamma = T)

# X correct, C incorrect
assess.dat.gamma(n = n, q = q, B = B, s2 = s2, x.means = x.means,
                 x.shape = 2, c.shape = 2,
                 specify.x.gamma = T, specify.c.gamma = F)

# X both Incorrect
assess.dat.gamma(n = n, q = q, B = B, s2 = s2, x.means = x.means,
                 x.shape = 2, c.shape = 2,
                 specify.x.gamma = F, specify.c.gamma = F)

# assess estimation of etas -----------------------------------------------

n.rep <- 100

assess.eta.est.gamma(n.rep = n.rep, n = n, q = q, B = B, s2 = s2,
                     x.means = x.means, x.shape = 2, c.shape = 2,
                     specify.x.gamma = T, specify.c.gamma = T)

