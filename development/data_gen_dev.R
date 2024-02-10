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
q <- 0.8              # censoring proportion
B <- c(1, 2)          # beta
s2 <- 1.1             # variance of Y|X,Z
x.mean <- 0.25        # mean of X

# generate data -----------------------------------------------------------

assess.dat(n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
           x.shape = 1, c.shape = 1)

assess.dat(n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
           x.shape = 1, c.shape = 4)

assess.dat(n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
           x.shape = 4, c.shape = 1)

assess.dat(n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
           x.shape = 4, c.shape = 4)
