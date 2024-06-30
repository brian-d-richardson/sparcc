###############################################################################
###############################################################################

# MLE B-spline Implementation

# Brian Richardson

# 2024-06-21

# Purpose: develop B-spline model for X|Z

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
library(statmod)
library(devtools)
library(ggplot2)
library(zipfR)
library(splines2)
library(nloptr)
library(tictoc)
load_all()

# define parameters -------------------------------------------------------

set.seed(3)
n <- 8000             # sample size
q <- 0.8              # censoring proportion
B <- c(1, 10, 2)      # beta
s2 <- 2               # variance of Y|X,Z
ls2 <- log(2)
x.thetas <- c(0.5, -0.5)   # theta parameters for X|Z
x.gamma <- 1          # gamma parameter for X|Z
c.gamma <- 2          # gamma parameter for C|Z
mx <- 40                   # number of nodes for X
my <- 5                    # number of nodes for Y

# generate data -----------------------------------------------------------

dat.list <- gen.data.beta(
  n = n, q = q, B = B, s2 = s2,
  x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma)

datf <- dat.list$datf          # full data
dat0 <- dat.list$dat0          # oracle data
dat <- dat.list$dat            # observed data
datcc <- dat.list$datcc        # complete case data
zs <- sort(unique(dat$Z))      # unique z values

# plot generated data
ggplot(datf,
       aes(x = X,
           y = Y,
           color = factor(Z))) +
  geom_point() +
  labs(color = "Z") +
  ggtitle("Oracle Data")

ggplot(datcc,
       aes(x = W,
           y = Y,
           color = factor(Z))) +
  geom_point() +
  labs(color = "Z") +
  ggtitle("Complete Case Data")

# estimate X|Z density using Splines --------------------------------------

# parameters for spline model
m.knots <- 5
deg <- 3

# fit spline model and extract estimated density
tic("fit spline model")
spline.res <- fit.spline(dat = dat, m.knots = m.knots, deg = deg)
toc()
eta1.hat <- spline.res$dens
alpha.hat <- spline.res$alpha
theta.hat <- spline.res$theta
logf <- spline.res$glogf
glogf <- spline.res$glogf
wts <- spline.res$wts
knots <- spline.res$knots

# plot results
datf %>%
  mutate(eta1.hat = eta1.hat(x = X, z = Z),
         Z = factor(Z)) %>%
  ggplot() +
  geom_line(aes(x = X,
                y = eta1.hat,
                color = Z)) +
  geom_histogram(aes(x = X,
                     y = after_stat(density),
                     fill = Z),
                 alpha = 0.5,
                 bins = 30) +
  facet_grid(~Z) +
  theme(legend.position = "none")

datcc %>%
  mutate(eta1.hat = eta1.hat(x = W, z = Z),
         Z = factor(Z)) %>%
  ggplot() +
  geom_line(aes(x = W,
                y = eta1.hat,
                color = Z)) +
  geom_histogram(aes(x = W,
                     y = after_stat(density),
                     fill = Z),
                 alpha = 0.5,
                 bins = 30) +
  facet_grid(~Z) +
  theme(legend.position = "none")

# create quadrature rules -------------------------------------------------

# X|Z quadrature
x.nds <- vapply(
  X = zs,
  FUN.VALUE = numeric(mx),
  FUN = function(z) seq(1E-6, 1-1E-6, length = mx))

x.wts <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mx),
  FUN = function(i) eta1.hat(x.nds[,i], zs[i]) /
    sum(eta1.hat(x.nds[,i], zs[i])))

# Y quadrature
gq <- gauss.quad(n = my, kind = "hermite")
y.nds <- gq$nodes
y.wts <- gq$weights

# fit outcome model -------------------------------------------------------

# complete case lm to get starting value
naive.lm <- lm(Y ~ W + Z, data = datcc)

# complete case
Bcc <- get.root(dat = dat, score = get.Scc,
                start = c(naive.lm$coef, log(var(naive.lm$resid))))
Bcc

# oracle
B0 <- get.root(dat = dat0, score = get.Scc,
               start = Bcc)
B0

# MLE
Bmle <- get.root(dat = dat, score = get.Sml, start = Bcc,
                 args = list(x.nds = x.nds, x.wts = x.wts))
Bmle

# compare estimates
rbind(c(B, ls2), B0, Bcc, Bmle)

# estimate variance of MLE ------------------------------------------------

# check derivative of log-likelihood at estimate
colSums(glogf(dat = dat, theta = theta.hat))

# MLE variance including nuisance estimation
tic("MLE variance estimation (outer product method)")
V.ml <- var.est.mle(
  dat = dat,
  get.S = function(dat, theta, args) {

    # split parameter into outcome and nuisance params
    theta.out <- head(theta, 4)
    theta.nuis <- matrix(tail(theta, -4), ncol = length(zs))

    # define estimated X|Z density
    eta1 <- function(x, z) {

      # spline basis for supplied x values
      bs.x <- splines::bs(
        x = x,
        knots = knots,
        Boundary.knots = c(0, 1),
        degree = deg,
        intercept = F)

      # transform theta parameter to alpha
      alpha.hat <- apply(theta.nuis, 2, function(x) softmax(x, wts = wts))

      # estimated density
      fxz <- numeric(length(x))
      for (zi in 1:length(zs)) {
        z.ind <- z == zs[zi]
        fxz[z.ind] <- bs.x[z.ind,] %*% alpha.hat[,zi]
      }
      return(fxz)
    }

    # create quadrature nodes
    x.wts <- vapply(
      X = 1:length(zs),
      FUN.VALUE = numeric(mx),
      FUN = function(i) eta1(x.nds[,i], zs[i]) / sum(eta1(x.nds[,i], zs[i])))

    args <- list(x.nds = x.nds, x.wts = x.wts)

    # stack estimating equations
    cbind(get.Sml(dat = dat, theta = head(theta, 4),
              args = args, return.sums = F),
          glogf(dat = dat, theta = matrix(theta.nuis, ncol = length(zs))))

  },
  theta = c(Bmle, theta.hat),
  args = list(),
  ridge.size = 1E-6,
  return.se = F
)
toc()

10 * sqrt(diag(V.ml)[1:4])

# sandwich variance excluding nuisance estimation
tic("Naive sandwich variance estimation")
V.ml2 <- var.est.sand(
  dat = dat,
  get.S = function(dat, theta, args, return.sums = F) {

  get.Sml(dat = dat, theta = theta,
          args = args, return.sums = return.sums)

  },
  theta = Bmle,
  args = list(),
  return.se = F
)
toc()

10 * sqrt(diag(V.ml2)[1:4])

# sandwich variance including nuisance estimation
tic("Sandwich variance estimation")
V.ml3 <- var.est.sand(
  dat = dat,
  get.S = function(dat, theta, args, return.sums = F) {

    # split parameter into outcome and nuisance params
    theta.out <- head(theta, 4)
    theta.nuis <- matrix(tail(theta, -4), ncol = length(zs))

    # define estimated X|Z density
    eta1 <- function(x, z) {

      # spline basis for supplied x values
      bs.x <- splines::bs(
        x = x,
        knots = knots,
        Boundary.knots = c(0, 1),
        degree = deg,
        intercept = F)

      # transform theta parameter to alpha
      alpha.hat <- apply(theta.nuis, 2, function(x) softmax(x, wts = wts))

      # estimated density
      fxz <- numeric(length(x))
      for (zi in 1:length(zs)) {
        z.ind <- z == zs[zi]
        fxz[z.ind] <- bs.x[z.ind,] %*% alpha.hat[,zi]
      }
      return(fxz)
    }

    # create quadrature nodes
    x.wts <- vapply(
      X = 1:length(zs),
      FUN.VALUE = numeric(mx),
      FUN = function(i) eta1(x.nds[,i], zs[i]) / sum(eta1(x.nds[,i], zs[i])))

    args <- list(x.nds = x.nds, x.wts = x.wts)

    # stack estimating equations
    S <- cbind(get.Sml(dat = dat, theta = head(theta, 4),
                  args = args, return.sums = F),
          glogf(dat = dat, theta = matrix(theta.nuis, ncol = length(zs))))

    if (return.sums) {
      return(colSums(S))
    } else {
      return(S)
    }

  },
  theta = c(Bmle, theta.hat),
  args = list(),
  return.se = F
)
toc()

10 * sqrt(diag(V.ml3)[1:4])



