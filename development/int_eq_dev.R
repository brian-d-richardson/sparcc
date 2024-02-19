###############################################################################
###############################################################################

# Develop integral equation solver

# Brian Richardson

# 2024-01-19

# Purpose: develop integral equation solver

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
setwd(dirname(getwd()))
library(devtools)
load_all()

# define parameters -------------------------------------------------------

n <- 1000                  # sample size
s2 <- 0.5                  # variance of Y|X,Z
q <- .9                    # censoring proportion
B <- c(1, 2)               # outcome model parameters
mx <- 21                   # number of nodes for X
mc <- 21                   # number of nodes for C
my <- 11                   # number of nodes for Y
x.rate <- 4                # rate parameter for gamma distribution of X
c.rate <- x.rate*q/(1-q)   # rate parameter for exp distribution of C

# X density
eta1 <- function(x) dexp(x, rate = x.rate)

# C density
eta2 <- function(c) dexp(c, rate = c.rate)

# mean function mu(X, B) = E(Y | X)
mu <- function(x, B) {
  B[1] + B[2]*x
}

# gradient of mu w.r.t. B
d.mu <- function(x, B) {
  cbind(1, x)
}

# Y density
fy <- function(y, x, B, s2) dnorm(x = y, mean = mu(x, B), sd = sqrt(s2))

# full data score vector
SF <- function(y, x, B, s2) {
  cbind((y - mu(x, B)) * d.mu(x, B),
        (y - mu(x, B)) ^ 2 - s2)
}

# generate data -----------------------------------------------------------

X <- rexp(n, x.rate)                            # censored covariate
C <- rexp(n, c.rate)                            # censoring time
W <- ifelse(X <= C, X, C)                       # observed covariate
Delta <- ifelse(X <= C, 1, 0)                   # uncensored indicator
Y <- rnorm(n, cbind(1, X) %*% B, sd = sqrt(s2)) # outcome
dat0 <- data.frame(Y, W = X, Delta = 1)         # oracle data
dat <- data.frame(Y, W, Delta)                  # observed data
datcc <- dat[Delta == 1,]                       # complete case data

# create quadrature rules -------------------------------------------------

# X quadrature
x.nds <- seq(10^-6, max(X), length = mx)
x.wts <- eta1(x.nds) / sum(eta1(x.nds))

# C quadrature
c.nds <- seq(10^-6, max(C), length = mc)
c.wts <- eta2(c.nds) / sum(eta2(c.nds))

# Y quadrature
gq <- gauss.quad(n = my, kind = "hermite")
y.nds <- gq$nodes
y.wts <- gq$weights

# solve integral equation -------------------------------------------------

# solve integral equation for a at node points
a <- get.a(B = B, s2 = s2, mu = mu, d.mu = d.mu, fy = fy, SF = SF,
           x.nds = x.nds, x.wts = x.wts, c.nds = c.nds, c.wts = c.wts,
           y.nds = y.nds, y.wts = y.wts)

# check that a has mean zero
colSums(a * x.wts)

# plot a
adat <- data.frame(cbind(x.nds, a, a*x.wts)) |>
  `colnames<-`(c("x", "a1", "a2", "a3", "a.eta1", "a.eta2", "a.eta3")) |>
  pivot_longer(cols = !x) |>
  mutate(comp = gsub("\\D", x = name, replacement = ""),
         func = gsub("\\d", x = name, replacement = ""))

ggplot(adat,
       aes(x = x,
           y = value,
           color = comp)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ func,
             scales = "free")
