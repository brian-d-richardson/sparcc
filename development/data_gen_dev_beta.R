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

# function to assess generated data ---------------------------------------

assess.dat.beta <- function(n, q, B, s2, x.thetas, x.gamma, c.gamma,
                            x.correct, c.correct) {

  ## generate data
  dat.list <- gen.data.beta(n = n, q = q, B = B, s2 = s2,
                            x.thetas = x.thetas,
                            x.gamma = x.gamma, c.gamma = c.gamma)

  datf <- dat.list$datf          # full data
  dat0 <- dat.list$dat0          # oracle data
  dat <- dat.list$dat            # observed data
  datcc <- dat.list$datcc        # complete case data

  ## compute empirical values by Z level
  q.hats <- datf %>% group_by(Z) %>% summarise(qhat = mean(X > C)) # cens prop
  x.bars <- datf %>% group_by(Z) %>% summarise(xbar = mean(X))     # mean X

  # estimate distribution of X|Z
  if (x.correct == T) {

    # estimate beta parameters at each level of Z
    x.params.hat <- t(vapply(
      X = 0:1,
      FUN.VALUE = numeric(2),
      FUN = function(z) {
        est <- dat %>%
          filter(Z == z) %>%
          mutate(left = W,
                 right = ifelse(Delta == 1, W, NA)) %>%
          dplyr::select(left, right) %>%
          fitdistrplus::fitdistcens(distr = "beta")
        return(est$estimate)
      }))

    # define estimated X|Z density
    eta1 <- function(x, z) {
      dbeta(x = x,
            shape1 = x.params.hat[z + 1, "shape1"],
            shape2 = x.params.hat[z + 1, "shape2"])
    }

  } else {

    # misspecify: estimate marginal beta distribution
    est <- dat %>%
      mutate(left = W,
             right = ifelse(Delta == 1, W, NA)) %>%
      dplyr::select(left, right) %>%
      fitdistrplus::fitdistcens(distr = "beta")
    x.params.hat <- est$estimate

    # define estimated X|Z density
    eta1 <- function(x, z) {
      dbeta(x = x,
            shape1 = x.params.hat["shape1"],
            shape2 = x.params.hat["shape2"])
    }
  }

  # estimate distribution of C|Z
  if (c.correct) {

    # estimate beta parameters at each level of Z
    c.params.hat <- t(vapply(
      X = 0:1,
      FUN.VALUE = numeric(2),
      FUN = function(z) {
        est <- dat %>%
          filter(Z == z) %>%
          mutate(left = W,
                 right = ifelse(Delta == 0, W, NA)) %>%
          dplyr::select(left, right) %>%
          fitdistrplus::fitdistcens(distr = "beta")
        return(est$estimate)
      }))

    # define estimated C|Z density
    eta2 <- function(c, z) {
      dbeta(x = c,
            shape1 = c.params.hat[z + 1, "shape1"],
            shape2 = c.params.hat[z + 1, "shape2"])
    }

  } else {

    # misspecify: estimate marginal beta distribution
    est <- dat %>%
      mutate(left = W,
             right = ifelse(Delta == 0, W, NA)) %>%
      dplyr::select(left, right) %>%
      fitdistrplus::fitdistcens(distr = "beta")
    c.params.hat <- est$estimate

    # define estimated X|Z density
    eta2 <- function(c, z) {
      dbeta(x = c,
            shape1 = c.params.hat["shape1"],
            shape2 = c.params.hat["shape2"])
    }
  }

  # grid to plot X density
  x.grid <- seq(0, 1, length = 100)

  # plot data
  datf %>%
    mutate(e1 = eta1(x = X, z = Z),
           e2 = eta2(c = C, z = Z)) %>%
    ggplot() +
    geom_histogram(aes(x = X,
                       y = after_stat(density)),
                   fill = "blue", alpha = 0.5, bins = n/50) +
    geom_histogram(aes(x = C,
                       y = after_stat(density)),
                   fill = "red", alpha = 0.5, bins = n/50) +
    geom_line(aes(x = X,
                  y = e1),
              color = "blue", linewidth = 1) +
    geom_line(aes(x = C,
                  y = e2),
              color = "red", linewidth = 1) +
    facet_wrap(~ Z,
               scales = "free", labeller = label_both) +
    labs(x = "X (blue) or C (red)",
         y = "Density") +
    ggtitle("Estimated vs Observed Distributions of X|Z and C|Z",
            subtitle = paste0("q.hats = ", paste(round(q.hats$qhat, 2), collapse = ", "), "; ",
                              "x.bars = ", paste(round(x.bars$xbar, 2), collapse = ", ")))

}

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

