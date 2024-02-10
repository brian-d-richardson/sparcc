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
install_github("brian-d-richardson/sparcc")
library(sparcc)
library(ggplot2)

# define parameters -------------------------------------------------------

n <- 10000            # sample size
q <- 0.4              # censoring proportion
B <- c(1, 2)          # beta
x.mean <- 0.25        # mean of X

# generate data -----------------------------------------------------------

assess.dat(n = n, q = q, B = B, x.mean = x.mean, x.shape = 1, c.shape = 1)
assess.dat(n = n, q = q, B = B, x.mean = x.mean, x.shape = 1, c.shape = 4)
assess.dat(n = n, q = q, B = B, x.mean = x.mean, x.shape = 4, c.shape = 1)
assess.dat(n = n, q = q, B = B, x.mean = x.mean, x.shape = 4, c.shape = 4)
