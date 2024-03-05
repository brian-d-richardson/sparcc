
# sparcc: semiparametric censored covariate estimation <img id="sparcc_hex" src="man/figures/sparcc_hex.png" align="right" width="125"/>

Brian D. Richardson

## Installation

Installation of `sparcc` from GitHub requires the
[`devtools`](https://www.r-project.org/nosvn/pandoc/devtools.html)
package and can be done with the following code:

``` r
# install the package
devtools::install_github(repo = "brian-d-richardson/sparcc", 
                         ref = "main")
```

``` r
# load the package
library(sparcc)

# other necessary packages
library(conf)
```

The `sparcc` package contains functions to analyze data with a random
right censored covariate using a semiparametric method. The methods
implemented are introduced in the paper, “Doubly robust estimation under
a random right censored covariate,” which is currently in progress.

## Example

Below is an example of three estimation procedures used on a data set
with a censored covariate. We first define parameters for data
generation and estimation.

``` r
# define parameters -------------------------------------------------------

set.seed(1)                # random number seed for reproducibility
B <- c(1, 2)               # outcome model parameters
s2 <- 1.1                  # Var(Y|X,Z)
q <- 0.6                   # censoring proportion
n <- 5000                  # sample size
x.mean <- 0.5              # mean of X
x.shape <- 1.1             # gamma shape parameter for X
c.shape <- 1.5             # gamma shape parameter for C
x.rate <- x.shape / x.mean # rate parameter for gamma distribution of X
c.rate <- get.c.rate(      # rate parameter for gamma distribution of C
  q = q,
  x.rate = x.rate,
  x.shape = x.shape,
  c.shape = c.shape)
mx <- 100                  # nodes in quadrature grid for X
mc <- 15                   # nodes in quadrature grid for C
my <- 3                    # nodes in quadrature grid for Y
specify.x.gamma <- F       # indicator for estimating X as gamma
specify.c.gamma <- T       # indicator for estimating C as gamma

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
```

We now generate data using the built in `gen.data` function. This
returns a list of data frames: (i) `datf`: The full data, including the
outcome `Y`, the covariate `X`, and the censoring time `C`. (ii) `dat`:
The observed data, including the outcome `Y`, the possibly censored
covariate `W`, and the censoring indicator `Delta`. (iii) `dat0`: The
oracle data, a version of the observed data where no observations are
censored (essentially setting `C` equal to infinity for all
observations) (iv) `datcc`: The complete case data, or the subset of
`dat` with `Delta == 1`.

``` r
dat.list <- gen.data(n = n, q = q, B = B, s2 = s2,
                     x.rate = x.rate, x.shape = x.shape,
                     c.rate = c.rate, c.shape = c.shape)
datf <- dat.list$datf          # full data
dat0 <- dat.list$dat0          # oracle data
dat <- dat.list$dat            # observed data
datcc <- dat.list$datcc        # complete case data
```

Using the observed data, we can estimated the nuisance distributions of
$X$ and $C$. If the prespecified parameter `specify.x.gamma` is `TRUE`,
then we (correctly) model $X$ with a gamma distribution. Else, we
(incorrectly) model $X$ with an exponential distribution. The same holds
for specifying a model for $C$ via `specify.c.gamma`.

In this example, we incorrectly model $X$ and correctly model $C$.

``` r
# estimate distribution of X
if (specify.x.gamma) {
  x.param.hat <- gammaMLE(yi = dat$W,                      # gamma parameters
                          si = dat$Delta,
                          scale = F)$estimate
  eta1 <- function(x) dgamma(x = x,                        # gamma density
                             shape = x.param.hat["shape"],
                             rate = x.param.hat["rate"])
} else {
  x.rate.hat <- mean(dat$Delta) / mean(dat$W)      # exponential rate parameter
  eta1 <- function(x) dexp(x, rate = x.rate.hat)   # exponential density
}

# estimate distribution of C
if (specify.c.gamma) {
  c.param.hat <- gammaMLE(yi = dat$W,                      # gamma parameters
                          si = 1 - dat$Delta,
                          scale = F)$estimate
  eta2 <- function(x) dgamma(x = x,                        # gamma density
                             shape = c.param.hat["shape"],
                             rate = c.param.hat["rate"])
} else {
  c.rate.hat <- mean(1 - dat$Delta) / mean(dat$W)  # exponential rate parameter
  eta2 <- function(x) dexp(x, rate = c.rate.hat)   # exponential density
}
```

Using the estimated densities, we then create quadrature rules for $X$
and $C$. Additionally, we create a Gauss-Hermite quadrature rule for
$Y$.

``` r
x.upper <- max(datf$X)  # oracle max
c.upper <- max(datf$C)

# X quadrature
x.nds <- seq(10^-5, x.upper, length = mx)
x.wts <- eta1(x.nds) / sum(eta1(x.nds))

# C quadrature
c.nds <- seq(10^-5, c.upper, length = mc)
c.wts <- eta2(c.nds) / sum(eta2(c.nds))

# Y quadrature
gq <- statmod::gauss.quad(n = my, kind = "hermite")
y.nds <- gq$nodes
y.wts <- gq$weights
```

We can now estimate the parameters of interest in the outcome model
using:

1)  oracle method (using the unobserved $X$ values)
2)  complete case
3)  maximum likelihood estimation
4)  semiparametric efficient estimation

Note that, since we misspecified $X$, the MLE is not a consistent
estimator. The semiparametric method on the other hand remains
consistent since the distribution for $C$ is correctly specified as
gamma.

``` r
# complete case lm to get starting value
naive.lm <- lm(Y ~ W, data = datcc)

# complete case
Bcc <- get.root(dat = dat, score = get.Scc,
                start = c(naive.lm$coef, log(var(naive.lm$resid))))

# oracle
B0 <- get.root(dat = dat0, score = get.Scc,
               start = Bcc)

# MLE
Bmle <- get.root(dat = dat, score = get.Sml, start = Bcc,
                 args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                             x.nds = x.nds, x.wts = x.wts))

# semiparametric efficient score
Beff <- get.root(dat = dat, score = get.Seff, start = Bcc,
                 args = list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                             eta1 = eta1,
                             x.nds = x.nds, x.wts = x.wts,
                             c.nds = c.nds, c.wts = c.wts,
                             y.nds = y.nds, y.wts = y.wts))
```

We then compare estimates. Note that the MLE estimate of the coefficient
for $X$ is the furthest from the truth out of the four estimators. This
makes sense since the MLE is inconsistent under a misspecified model for
$X$. Note also that the complete case is further from the truth than the
oracle and the semiparametric.

``` r
# compare estimates
rbind(c(B, log(s2)), B0, Bcc, Bmle, Beff)
```

    ##      (Intercept)        W           
    ##        1.0000000 2.000000 0.09531018
    ## B0     0.9832562 2.033754 0.07131662
    ## Bcc    0.9541664 2.042182 0.09957880
    ## Bmle   1.0096010 1.888328 0.06340745
    ## Beff   0.9601168 2.005944 0.09961160
