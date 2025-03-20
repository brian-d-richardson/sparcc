###############################################################################
###############################################################################

# Wrapper development

# Anonymous

# 2025-03-19

# Purpose: develop wrapper function for sparcc estimator

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
#library(statmod)
library(devtools)
#library(ggplot2)
#library(zipfR)
#library(splines2)
#library(nloptr)
#library(tictoc)
load_all()


# define parameters -------------------------------------------------------

seed <- 1                  # random number seed
n <- 800                   # sample size
q <- 0.4                   # censoring proportion

B <- c(1, 10, 2)           # outcome model parameters
s2 <- 1

x.thetas <- 0.5 * c(-1, 1) # nuisance model parameters
x.gamma <- 1
c.gamma <- 2

# generate data -----------------------------------------------------------

set.seed(seed)
dat.list <- gen.data.beta(
  n = n, q = q, B = B, s2 = s2,
  x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma)

data <- dat.list$dat            # observed data

# assess data -------------------------------------------------------------

# plot full data
ggplot(datf,
       aes(x = X,
           y = Y,
           color = factor(Z))) +
  geom_point() +
  labs(color = "Z") +
  ggtitle("Oracle Data")

# plot complete cases only
ggplot(datcc,
       aes(x = W,
           y = Y,
           color = factor(Z))) +
  geom_point() +
  labs(color = "Z") +
  ggtitle("Complete Case Data")

# fit sparcc estimator with parametric models -----------------------------

sparcc.param = sparcc(
  data = data,
  nuisance.models = "parametric",
  distr.x = "beta",
  distr.c = "beta"
)

# fit sparcc estimator with nonparametric models --------------------------

sparcc.nonpar = sparcc(
  data = data,
  nuisance.models = "nonparametric"
)





