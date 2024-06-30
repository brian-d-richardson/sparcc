###############################################################################
###############################################################################

# MLE Variance Estimation

# Brian Richardson

# 2024-06-18

# Purpose: develop variance estimator for MLE when X|Z ~ Beta

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
library(statmod)
library(devtools)
library(fitdistrplus)
library(ggplot2)
library(zipfR)
load_all()

# define parameters -------------------------------------------------------

set.seed(2)
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

# estimate nuisance distributions -----------------------------------------

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

# create quadrature rules -------------------------------------------------

# X|Z quadrature
x.nds <- vapply(
  X = zs,
  FUN.VALUE = numeric(mx),
  FUN = function(z) seq(1E-6, 1-1E-6, length = mx))

x.wts <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mx),
  FUN = function(i) eta1(x.nds[,i], zs[i]) / sum(eta1(x.nds[,i], zs[i])))

# Y quadrature
gq <- gauss.quad(n = my, kind = "hermite")
y.nds <- gq$nodes
y.wts <- gq$weights

# estimate beta -----------------------------------------------------------

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

# find X|Z estimating function --------------------------------------------

# check that derivative of X|Z log-likelihood is close to zero
dlf <- d.log.fxz(dat = dat, theta = c(x.params.hat), return.sums = F)
colSums(dlf)

# estimate variance of X|Z parameters
x.params.se <- var.est.sand(dat = dat,
                            get.S = d.log.fxz,
                            theta = c(x.params.hat),
                            return.se = T)

# compare to std error from fitdistrplus function
x.params.sd <- t(vapply(
  X = 0:1,
  FUN.VALUE = numeric(2),
  FUN = function(z) {
    est <- dat %>%
      filter(Z == z) %>%
      mutate(left = W,
             right = ifelse(Delta == 1, W, NA)) %>%
      dplyr::select(left, right) %>%
      fitdistrplus::fitdistcens(distr = "beta")
    return(est$sd)
  }))

# expect these to be very similar
x.params.se
c(x.params.sd)

# estimate variance of MLE parameters -------------------------------------

# sandwich variance ignoring nuisance estimation
V.ml.naive <- var.est.sand(
  dat = dat,
  get.S = get.Sml,
  theta = Bmle,
  args = list(x.nds = x.nds, x.wts = x.wts),
  return.se = F
)

# sandwich variance including nuisance estimation
V.ml <- var.est.sand(
  dat = dat,
  get.S = function(dat, theta, args, return.sums = F) {

    alpha <- tail(theta, -4)

    # define estimated X|Z density
    eta1 <- function(x, z) {
      dbeta(x = x,
            shape1 = alpha[z + 1],
            shape2 = alpha[z + 3])
    }

    # create quadrature nodes
    x.wts <- vapply(
      X = 1:length(zs),
      FUN.VALUE = numeric(mx),
      FUN = function(i) eta1(x.nds[,i], zs[i]) / sum(eta1(x.nds[,i], zs[i])))

    args <- list(x.nds = x.nds, x.wts = x.wts)

    # stack estimating equations
    cbind(get.Sml(dat = dat, theta = head(theta, 4),
                  args = args, return.sums = return.sums),
          d.log.fxz(dat = dat, theta = alpha,
                    args = args, return.sums = return.sums))
  },
  theta = c(Bmle, x.params.hat),
  args = list(x.nds = x.nds, x.wts = x.wts),
  return.se = F
)

# compare variance estimators (expect the naive version to be smaller)
sqrt(diag(V.ml.naive))
sqrt(diag(V.ml)[1:4])






