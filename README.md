
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
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

The `sparcc` package contains functions to analyze data with a random
right censored covariate using a semiparametric method. The methods
implemented are introduced in the paper, “Doubly robust estimation under
a random right censored covariate,” which is currently in progress.

The code implemented in this package is specific to the scenario where
$Y|X,Z$ has a normal distribution with mean
$\textrm{E}(Y|X,Z)=\beta_0+\beta_1X+\beta_2Z$, where $X$ is a censored
covariate and $Z$ is an uncensored covariate.

## Example

Below is an example of three estimation procedures used on a data set
with a censored covariate. We first define parameters for data
generation and estimation.

``` r
# define parameters -------------------------------------------------------

set.seed(123)                 # random number seed for reproducibility
n <- 8000                     # sample size
q <- 0.8                      # censoring proportion
B <- c(1, 10, 2)              # outcome model parameters
s2 <- 1                       # Var(Y|X,Z)
x.thetas <- 0.5 * c(-1, 1)    # parameters governing X|Z and C|Z
x.gamma <- 1
c.gamma <- 2
mx <- 40                      # number of nodes in quadrature grid for X|Z
mc <- 40                      # number of nodes in quadrature grid for C|Z
my <- 5                       # number of nodes in quadrature grid for Y|X,Z
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

For this example, $X|Z$ and $C|Z$ follow beta distributions leading to
the desired censoring proportion $q=\textrm{P}(X>C)$.

``` r
dat.list <- gen.data.beta(n = n, q = q, B = B, s2 = s2,
                          x.thetas = x.thetas, x.gamma = x.gamma, c.gamma)
datf <- dat.list$datf          # full data
dat0 <- dat.list$dat0          # oracle data
dat <- dat.list$dat            # observed data
datcc <- dat.list$datcc        # complete case data
zs <- sort(unique(dat$Z))      # unique z values
```

Using the observed data, we can estimate the nuisance distributions of
$X|Z$ and $C|Z$. To illustrate the consequences of model
misspecification, we incorrectly model $X|Z$ with marginal beta
distribution, while it truly follows a conditional beta distribution. We
correctly specify $C|Z$ as conditional beta.

``` r
## estimated parameters for X|Z
x.fit <- dat %>%
  mutate(left = W,
           right = ifelse(Delta == 1, W, NA)) %>%
  dplyr::select(left, right) %>%
  fitdistrplus::fitdistcens(distr = "beta")
x.params.hat <- x.fit$estimate

## estimated density of X|Z
eta1 <- function(x, z) {
  dbeta(x = x,
        shape1 = x.params.hat["shape1"],
        shape2 = x.params.hat["shape2"])
}

## estimated parameters for C|Z
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

## estimated density of C|Z
eta2 <- function(c, z) {
  dbeta(x = c,
        shape1 = c.params.hat[z + 1, "shape1"],
        shape2 = c.params.hat[z + 1, "shape2"])
}
```

Using the estimated densities, we then create quadrature rules for $X|Z$
and $C|Z$. Additionally, we create a Gauss-Hermite quadrature rule for
$Y$.

``` r
## X|Z quadrature nodes
x.nds <- vapply(
  X = zs,
  FUN.VALUE = numeric(mx),
  FUN = function(z) seq(1E-6, 1-1E-6, length = mx))

## X|Z quadrature weights
x.wts <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mx),
  FUN = function(i) eta1(x.nds[,i], zs[i]) /
    sum(eta1(x.nds[,i], zs[i])))

## C|Z quadrature nodes
c.nds <- vapply(
  X = zs,
  FUN.VALUE = numeric(mc),
  FUN = function(z) seq(1E-6, 1-1E-6, length = mc))

## X|Z quadrature weights
c.wts <- vapply(
  X = 1:length(zs),
  FUN.VALUE = numeric(mc),
  FUN = function(i) eta2(c.nds[,i], zs[i]) /
    sum(eta2(c.nds[,i], zs[i])))

## Y|X,Z quadrature
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

Note that, since we misspecified $X|Z$, the MLE is not a consistent
estimator. The semiparametric method on the other hand remains
consistent since the distribution for $C|Z$ is correctly specified as
gamma.

``` r
## complete case linear model to get starting value
naive.lm <- lm(Y ~ W + Z, data = datcc)

## complete case
B.cc <- get.root(dat = dat, score = get.Scc,
                 start = c(naive.lm$coef, log(var(naive.lm$resid))))

round(B.cc, 2)
```

    ##    B1    B2    B3    B4 
    ##  0.94 10.03  2.08  0.03

``` r
## oracle
B.or <- get.root(dat = dat0, score = get.Scc, start = B.cc)

round(B.or, 2)
```

    ##    B1    B2    B3    B4 
    ##  0.91 10.09  2.07  0.00

The MLE requires the additional arguments `x.nds` and `x.wts`,
corresponding to the nodes and weights in the quadrature rule for $X|Z$.
These are passed to the `get.root` function as named items in the `args`
list.

``` r
## MLE
mle.args <- list(x.nds = x.nds, x.wts = x.wts)
B.mle <- get.root(dat = dat, score = get.Sml, start = B.cc,
                  args = mle.args)

round(B.mle, 2)
```

    ##   B1   B2   B3   B4 
    ## 1.94 8.24 0.93 0.22

The MLE requires the additional arguments `x.nds`, `x.wts`, `c.nds`, and
`c.wts`, corresponding to the nodes and weights in the quadrature rules
for $X|Z$ and $C|Z$.

``` r
## semiparametric efficient estimator
sp.args <- list(x.nds = x.nds, x.wts = x.wts,
                c.nds = c.nds, c.wts = c.wts,
                y.nds = y.nds, y.wts = y.wts)
B.sp <- get.root(dat = dat, score = get.Seff, start = B.cc,
                 args = sp.args)

round(B.sp, 2)
```

    ##   B1   B2   B3   B4 
    ## 0.96 9.95 2.07 0.03

We then compare estimates. Note that the MLE estimates are far from the
truth,. which sense since the MLE is inconsistent under a misspecified
model for $X|Z$.

``` r
# compare estimates
round(rbind(c(B, log(s2)), B.or, B.cc, B.mle, B.sp), 2)
```

    ##         B1    B2   B3   B4
    ##       1.00 10.00 2.00 0.00
    ## B.or  0.91 10.09 2.07 0.00
    ## B.cc  0.94 10.03 2.08 0.03
    ## B.mle 1.94  8.24 0.93 0.22
    ## B.sp  0.96  9.95 2.07 0.03

Finally, we can compute standard errors for the different estimators
based on influence function variance estimators. For the MLE, this type
of variance estimation is not valid and so is not included.

``` r
# complete case
SE.cc <- var.est(dat = datcc, theta = B.cc, args = list(),
                 get.S = get.Scc, return.se = T)

## oracle
SE.or <- var.est(dat = dat0, theta = B.or, args = list(),
                 get.S = get.Scc, return.se = T)

## semiparametric efficient
SE.sp <- var.est(dat = dat, theta = B.sp,
                 args = sp.args, get.S = get.Seff, return.se = T)

round(rbind(SE.or, SE.cc, SE.sp), 3)
```

    ##        [,1]  [,2]  [,3]  [,4]
    ## SE.or 0.036 0.052 0.026 0.016
    ## SE.cc 0.080 0.182 0.066 0.037
    ## SE.sp 0.079 0.181 0.066 0.037
