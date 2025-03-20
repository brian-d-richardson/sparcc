###############################################################################
###############################################################################

# Develop integral equation solver

# Anonymous

# 2024-01-19

# Purpose: develop integral equation solver

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
library(devtools)
library(microbenchmark)
load_all()

# define parameters -------------------------------------------------------

seed <- 1
n <- 8000
q <- 0.8
B <- c(1, 10, 2)
s2 <- 1
x.thetas <- 0.5 * c(-1, 1)
x.gamma <- 1
c.gamma <- 2
mx <- 40
mc <- 40
my <- 5

# generate data -----------------------------------------------------------

set.seed(seed)
dat.list <- gen.data.beta(
  n = n, q = q, B = B, s2 = s2,
  x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma)

datf <- dat.list$datf          # full data
dat0 <- dat.list$dat0          # oracle data
dat <- dat.list$dat            # observed data
datcc <- dat.list$datcc        # complete case data
zs <- sort(unique(dat$Z))      # unique z values

# estimate nuisance distributions -----------------------------------------

## X|Z quadrature nodes
x.nds <- vapply(
  X = zs,
  FUN.VALUE = numeric(mx),
  FUN = function(z) seq(1E-6, 1-1E-6, length = mx))

## true distribution of X|Z
eta1 <- function(x, z) {
  dbeta(x = x,
        shape1 = 1 + x.gamma - x.thetas[z + 1],
        shape2 = 1 + x.gamma + x.thetas[z + 1])
}

## quadrature weights for X|Z
x.wts <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mx),
  FUN = function(i) eta1(x.nds[,i], zs[i]) /
    sum(eta1(x.nds[,i], zs[i])))

## theta parameters for beta distribution of C|Z
c.thetas <- vapply(
  X = x.thetas,
  FUN.VALUE = 0,
  FUN = function(xt)
    get.c.param.beta(q = q, x.theta = xt,
                     x.gamma = x.gamma, c.gamma = c.gamma))

## C|Z quadrature nodes
c.nds <- vapply(
  X = zs,
  FUN.VALUE = numeric(mc),
  FUN = function(z) seq(1E-6, 1-1E-6, length = mc))

## true distribution of C|Z
eta2 <- function(c, z) {
  dbeta(x = c,
        shape1 = 1 + c.gamma - c.thetas[z + 1],
        shape2 = 1 + c.gamma + c.thetas[z + 1])
}

c.wts <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mc),
  FUN = function(i) eta2(c.nds[,i], zs[i]) /
    sum(eta2(c.nds[,i], zs[i])))

## Y|X,Z quadrature nodes
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
  y.nds = y.nds, y.wts = y.wts,
  xz.interaction = F)

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
  tidyr::pivot_longer(cols = !c(x, z)) |>
  dplyr::mutate(comp = gsub("\\D", x = name, replacement = ""),
         func = gsub("\\d", x = name, replacement = ""),
         Function = ifelse(func == "a", "a", "a * eta1"))

ggplot(adat,
       aes(x = x,
           y = value,
           color = comp)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  facet_grid(Function ~ z,
             scales = "free",
             labeller = label_both) +
  labs(x = "X",
       y = "Function Value",
       color = "Component")

# Evaluate Computation Time -----------------------------------------------

solve.a <- function() {
  get.a(
    B = B, s2 = s2, mu = mu, d.mu = d.mu,
    fy = fy, SF = SF, zs = zs,
    x.nds = x.nds, x.wts = x.wts,
    c.nds = c.nds, c.wts = c.wts,
    y.nds = y.nds, y.wts = y.wts,
    xz.interaction = F)
}
microbenchmark(solve.a(), times = 10)
