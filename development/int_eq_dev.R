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
#setwd(dirname(getwd()))
library(devtools)
load_all()

# define parameters -------------------------------------------------------

n <- 1000                  # sample size
s2 <- 0.81                 # variance of Y|X,Z
q <- .8                    # censoring proportion
B <- c(1, 2, -1)           # outcome model parameters
x.means <- c(0.5, 1)       # mean of X at levels of Z
x.shape <- 2               # shape parameter for X
c.shape <- 2               # shape parameter for C
mx <- 30                   # number of nodes for X
mc <- 15                   # number of nodes for C
my <- 5                    # number of nodes for Y

x.rates <- x.shape / x.means  # rate parameters for gamma distribution of X|Z
c.rates <- vapply(
  X = x.rates,                # rate parameters for gamma distribution of C|Z
  FUN.VALUE = 0,
  FUN = function(xr)
    get.c.rate(q = q, x.rate = xr,
               x.shape = x.shape,
               c.shape = c.shape))

# X|Z density
eta1 <- function(x, z) {
  dgamma(x = x,
         shape = x.shape,
         rate = x.rates[z + 1])
}

# C|Z density
eta2 <- function(c, z) {
  dgamma(x = c,
         shape = c.shape,
         rate = c.rates[z + 1])
}

# mean function mu(X, Z, B) = E(Y | X, Z)
mu <- function(x, z, B) {
  B[1] + B[2]*x + B[3]*z
}

# gradient of mu w.r.t. B
d.mu <- function(x, z, B) {
  cbind(1, x, z)
}

# Y|X, Z density
fy <- function(y, x, z, B, s2) dnorm(x = y, mean = mu(x, z, B), sd = sqrt(s2))

# full data score vector
SF <- function(y, x, z, B, s2) {
  cbind((y - mu(x, z, B)) * d.mu(x, z, B),
        (y - mu(x, z, B)) ^ 2 - s2)
}

# generate data -----------------------------------------------------------

dat.list <- gen.data(n = n, q = q, B = B, s2 = s2,
                     x.means = x.means, x.shape = x.shape, c.shape = c.shape)

datf <- dat.list$datf          # full data
dat0 <- dat.list$dat0          # oracle data
dat <- dat.list$dat            # observed data
datcc <- dat.list$datcc        # complete case data
zs <- sort(unique(dat$Z))      # unique z values

# create quadrature rules -------------------------------------------------

# X|Z quadrature
x.nds <- vapply(
  X = zs,
  FUN.VALUE = numeric(mx),
  FUN = function(z) seq(10^-6, max(datf$X[datf$Z == z]), length = mx))

x.wts <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mx),
  FUN = function(i) eta1(x.nds[,i], zs[i]) / sum(eta1(x.nds[,i], zs[i])))

# C|Z quadrature
c.nds <- vapply(
  X = zs,
  FUN.VALUE = numeric(mc),
  FUN = function(z) seq(10^-6, max(datf$C[datf$Z == z]), length = mc))

c.wts <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mc),
  FUN = function(i) eta2(c.nds[,i], zs[i]) / sum(eta2(c.nds[,i], zs[i])))

# Y quadrature
gq <- gauss.quad(n = my, kind = "hermite")
y.nds <- gq$nodes
y.wts <- gq$weights

# solve integral equation -------------------------------------------------

# solve integral equation for a at node points, for each unique Z value
a <- get.a(
  B = B, s2 = s2, mu = mu, d.mu = d.mu,
  fy = fy, SF = SF, zs = zs,
  x.nds = x.nds, x.wts = x.wts,
  c.nds = c.nds, c.wts = c.wts,
  y.nds = y.nds, y.wts = y.wts)

# check that a has mean zero
colMeans(a[[1]] * x.wts[,1])
colMeans(a[[2]] * x.wts[,2])

# plot a
adat <- data.frame(rbind(
  cbind(zs[1], x.nds[,1], a[[1]], a[[1]]*x.wts[,1]),
  cbind(zs[2], x.nds[,2], a[[2]], a[[2]]*x.wts[,2]))) |>
  `colnames<-`(c("z", "x",
                 "a1", "a2", "a3", "a4",
                 "a.eta1", "a.eta2", "a.eta3", "a.eta4")) |>
  pivot_longer(cols = !c(x, z)) |>
  mutate(comp = gsub("\\D", x = name, replacement = ""),
         func = gsub("\\d", x = name, replacement = ""))

ggplot(adat,
       aes(x = x,
           y = value,
           color = comp)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  facet_wrap(z ~ func,
             scales = "free",
             labeller = label_both)
