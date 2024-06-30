#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Fit spline model with censored outcome data
#'
#' @param dat data frame, the observed data
#' @param m.knots a positive integer, the number of knots for the spline model
#' @param deg a positive integer, the degree of the spline basis functions
#'
#' @return a list including
#' \itemize{
#' \item{`alpha`: a matrix, unconstrained spline model parameters}
#' \item{`theta`: a matrix, unconstrained spline model parameters}
#' \item{`dens`: a function, the estimated density}
#' \item{`logf`: a function, the log-likelihood}
#' \item{`glogf`: a function, gradient of the log-likelihood}
#' \item{`deg`: a positive integer, degree of the spline basis functions}
#' \item{`knots`: a numeric vector, knots used in the spline model}
#' \item{`wts`: a numeric vector, weights used in the theta-to-alpha transformation}
#' }
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.spline <- function(dat, m.knots, deg) {

  # unique z values
  zs <- sort(unique(dat$Z))

  # B-spline knots
  knots <- quantile(c(0, 1, dat$W), probs = seq(0, 1, length = m.knots))

  # B-spline basis with grid for integration
  xx <- seq(0, 1, length = 100)
  bs.grid <- splines::bs(
    x = xx,
    knots = knots,
    Boundary.knots = c(0, 1),
    degree = deg,
    intercept = F)
  #bs.grid <- bs.grid[, -c(ncol(bs.grid) - 1, ncol(bs.grid))]
  wts <- colMeans(bs.grid)

  # negative log density of W, Delta
  log.fwd.spline <- function(dat, theta) {

    alpha <- softmax(theta, wts = wts)
    delta <- dat$Delta
    w <- dat$W
    n <- nrow(dat)
    subint.width <- 1/length(xx)

    # B-spline basis with observed W
    bs.obs <- splines::bs(
      x = dat$W,
      knots = knots,
      Boundary.knots = c(0, 1),
      degree = deg,
      intercept = F)
    #bs.obs <- bs.obs[, -c(ncol(bs.obs) - 1, ncol(bs.obs))]

    f <- numeric(n)

    for (ii in 1:n) {
      if (delta[ii] == 1) {
        f[ii] <- bs.obs[ii,] %*% alpha
      } else {
        f[ii] <- subint.width * sum(bs.grid[xx > w[ii],] %*% alpha)
      }
    }
    return(-log(f))
  }

  # negative log-density of W, Delta | Z
  log.fwdz.spline <- function(dat, theta) {

    ll <- numeric(nrow(dat))
    for (zi in 1:length(zs)) {

      # subset to observations with given Z value
      z.ind <- dat$Z == zs[zi]
      datz <- dat[z.ind,]

      # compute loglik value
      ll[z.ind] <- log.fwd.spline(dat = datz, theta = theta[,zi])
    }
    return(ll)
  }

  # derivative of negative log density of W, Delta w.r.t. theta
  g.log.fwd.spline <- function(dat, theta) {

    alpha <- softmax(theta, wts = wts)

    delta <- dat$Delta
    w <- dat$W
    n <- nrow(dat)
    subint.width <- 1/length(xx)

    # B-spline basis with observed W
    bs.obs <- splines::bs(
      x = dat$W,
      knots = knots,
      Boundary.knots = c(0, 1),
      degree = deg,
      intercept = F)
    #bs.obs <- bs.obs[, -c(ncol(bs.obs) - 1, ncol(bs.obs))]

    # initiate f: density
    f <- numeric(n)

    # initiate g: derivative of f w.r.t. alpha
    g <- matrix(nrow = n, ncol = length(alpha))

    for (ii in 1:n) {
      if (delta[ii] == 1) {
        f[ii] <- bs.obs[ii,] %*% alpha
        g[ii,] <- bs.obs[ii,]
      } else {
        f[ii] <- subint.width * sum(bs.grid[xx > w[ii],] %*% alpha)
        if (sum(xx > w[ii]) == 1) {
          g[ii,] <- subint.width * bs.grid[xx > w[ii],]
        } else {
          g[ii,] <- subint.width * colSums(bs.grid[xx > w[ii],])
        }

      }
    }

    # apply chain rule to get derivative wrt theta
    gg <- (-g / f) %*% g.softmax(theta, wts = wts)

    return(gg)
  }

  # negative derivative of log-density of W, Delta | Z
  g.log.fwdz.spline <- function(dat, theta) {

    # initiate derivative matrix
    gll <- matrix(0, nrow = nrow(dat), ncol = length(c(theta)))

    for (zi in 1:length(zs)) {

      # subset to observations with given Z value
      z.ind <- dat$Z == zs[zi]
      z.yind <- 1:nrow(theta) + (zi - 1) * nrow(theta)
      datz <- dat[z.ind,]

      # compute loglik value
      gll[z.ind, z.yind] <-
        g.log.fwd.spline(dat = datz, theta = theta[,zi])

    }
    return(gll)
  }

  # initial parameters
  theta0 <- matrix(0, nrow = length(wts) - 1, ncol = length(zs))

  # estimate B-spline coefficients
  opt.res <- optim(
    par = c(theta0),
    fn = function(x) sum(log.fwdz.spline(
      dat = dat,
      theta = matrix(x, ncol = length(zs)))),
    gr = function(x) colSums(g.log.fwdz.spline(
      dat = dat,
      theta = matrix(x, ncol = length(zs)))),
    method = "BFGS", control = list(maxit = 1000)
  )

  theta.hat <- matrix(opt.res$par, ncol = length(zs))
  alpha.hat <- apply(theta.hat, 2,
                     function(t) softmax(x = t, wts = wts))

  # estimated density of X
  get.fxz.hat <- function(x, z) {

    # spline basis for supplied x values
    bs.x <- splines::bs(
      x = x,
      knots = knots,
      Boundary.knots = c(0, 1),
      degree = deg,
      intercept = F)
    #bs.x <- bs.x[, -c(ncol(bs.x) - 1, ncol(bs.x))]

    # estimated density
    fxz <- numeric(length(x))
    for (zi in 1:length(zs)) {
      z.ind <- z == zs[zi]
      fxz[z.ind] <- bs.x[z.ind,] %*% alpha.hat[,zi]
    }
    return(fxz)
  }

  # return parameters and density
  return(list(alpha = alpha.hat,
              theta = theta.hat,
              dens = get.fxz.hat,
              logf = log.fwdz.spline,
              glogf = g.log.fwdz.spline,
              deg = deg,
              knots = knots,
              wts = wts))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Softmax transformation
#'
#' @param x a numeric vector, values to be transformed
#' @param wts a numeric vector, weights
#'
#' @return a numeric vector
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
softmax <- function(x, wts) {
  c(exp(x), 1) / sum(wts * c(exp(x), 1))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Inverse softmax transformation
#'
#' @inheritParams softmax
#'
#' @return a numeric vector
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inv.softmax <- function(x, wts) {

  log(head(x, -1) / tail(x, 1))

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Derivative of softmax transformation
#'
#' @inheritParams softmax
#'
#' @return a numeric matrix
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g.softmax <- function(x, wts) {

  m <- length(x)
  denom <- sum(wts * c(exp(x), 1))

  # initialize derivative matrix
  g <- matrix(0, nrow = m + 1, ncol = m)

  # Compute the Jacobian matrix
  for (i in 1:m) {
    for (k in 1:m) {
      if (i == k) {
        g[i, k] <- (exp(x[k]) * (denom - wts[k] * exp(x[k]))) / denom ^ 2
      } else {
        g[i, k] <- (-wts[k] * exp(x[i] + x[k])) / denom ^ 2
      }
    }
  }

  # Last row for the constant term
  g[m + 1, ] <- -wts[k] * exp(x[k]) / denom ^ 2

  return(g)
}
