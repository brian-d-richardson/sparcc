###############################################################################
###############################################################################

# Develop Data Generation when X and C have beta distributions

# Brian Richardson

# 2024-01-22

# Purpose: develop methods to generate data

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
library(ggplot2)
library(dplyr)
library(devtools)
library(fitdistrplus)
library(dplyr)
load_all()

# define parameters -------------------------------------------------------

n <- 8000             # sample size
q <- 0.8              # censoring proportion
B <- c(1, 10, 2)      # beta
s2 <- 1               # variance of Y|X,Z
x.thetas <- c(0.5, -0.5)   # theta parameters for X|Z
x.gamma <- 1          # gamma parameter for X|Z
c.gamma <- 2          # gamma parameter for C|Z

# generate data -----------------------------------------------------------

# both correct
assess.dat.beta(n = n, q = q, B = B, s2 = s2,
                x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma,
                x.correct = T, c.correct = T)

# X incorrect, C correct
assess.dat.beta(n = n, q = q, B = B, s2 = s2,
                x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma,
                x.correct = F, c.correct = T)

# X correct, C incorrect
assess.dat.beta(n = n, q = q, B = B, s2 = s2,
                x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma,
                x.correct = T, c.correct = F)

# X both Incorrect
assess.dat.beta(n = n, q = q, B = B, s2 = s2,
                x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma,
                x.correct = F, c.correct = F)
