#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' compute d vector in integral equation
#'
#' @inheritParams get.b
#'
#' @return d matrix
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.d <- function(x.nds, x.wts, c.nds, c.wts) {

  mx <- length(x.nds)
  d <- numeric(mx)
  for (i in 1:mx) {
    d[i] <- sum(c.wts[x.nds[i] <= c.nds])
  }
  d <- d / x.wts
  return(d)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' compute M matrix in integral equation
#'
#' @inheritParams get.b
#'
#' @return M matrix
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.M <- function(B, s2, mu, fy, z,
                  x.nds, x.wts, c.nds, c.wts, y.nds, y.wts) {

  mx <- length(x.nds)
  mc <- length(c.nds)
  my <- length(y.nds)
  M <- matrix(0, nrow = mx, ncol = mx)
  for (i in 1:mx) {
    for (j in i:mx) {
      Mij <- 0
      mubar <- 0.5*(mu(x.nds[i], z, B) + mu(x.nds[j], z, B))
      for (h in 1:my) {
        for (l in 1:mc) {
          if (c.nds[l] < x.nds[i]) {
            Mij <- Mij + y.wts[h] * c.wts[l] /
            sum(fy(y = sqrt(s2)*y.nds[h] + mubar,
                   x = x.nds[x.nds > c.nds[l]],
                   z = z, B = B, s2 = s2) *
                x.wts[x.nds > c.nds[l]])
          }
        }
      }
      Mij <- Mij * exp(-0.25*(B[2]*(x.nds[i] - x.nds[j]))^2)
      M[i,j] <- M[j,i] <- Mij
    }
  }
  M <- M / (2*pi*sqrt(s2))
  return(M)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' compute b matrix in integral equation
#'
#' @param B a numeric vector, parameters in the outcome model
#' @param s2 a positive number, variance in the outcome model
#' @param mu a function, mean of outcome Y given covariates X, Z
#' @param fy a function, conditional density of outcome Y given covariates X, Z
#' @param SF a function, full data score
#' @param x.nds a numeric vector, nodes for quadrature rule for X
#' @param x.wts a numeric vector, weights for quadrature rule for X
#' @param c.nds a numeric vector, nodes for quadrature rule for C
#' @param c.wts a numeric vector, weights for quadrature rule for C
#' @param y.nds a numeric vector, nodes for quadrature rule for Y
#' @param y.wts a numeric vector, weights for quadrature rule for Y
#'
#' @return b matrix
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.b <- function(B, s2, mu, d.mu, fy, SF, z,
                  x.nds, x.wts, c.nds, c.wts, y.nds, y.wts) {

  mx <- length(x.nds)
  mc <- length(c.nds)
  my <- length(y.nds)
  b <- matrix(0, nrow = length(B) + 1, ncol = mx)
  for (i in 1:mx) {
    bi <- 0
    mui <- mu(x.nds[i], z, B)
    for (h in 1:my) {
      yih <- sqrt(2*s2)*y.nds[h] + mui
      for (l in 1:mc) {
        if (c.nds[l] < x.nds[i]) {
          xgc <- x.nds > c.nds[l]
          bi <- bi + y.wts[h] * c.wts[l] *
            colSums(fy(y = yih, x = x.nds[xgc], z = z, B = B, s2 = s2) *
                    x.wts[xgc] *
                    SF(y = yih, x = x.nds[xgc], z = z, B = B, ls2 = log(s2))) /
            sum(fy(y = yih, x = x.nds[xgc], z = z, B = B, s2 = s2) *
                x.wts[xgc])
        }
      }
    }
    b[,i] <- bi
  }
  b <- b / sqrt(pi)
  return(b)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' compute a function at nodes for a single level of z
#'
#' @inheritParams get.b
#'
#' @return approximated a values at nodes
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.a.z <- function(B, s2, mu, d.mu, fy, SF, z,
                    x.nds, x.wts, c.nds, c.wts, y.wts, y.nds) {

  d <- get.d(x.nds = x.nds, x.wts = x.wts,
             c.nds = c.nds, c.wts = c.wts)
  M <- get.M(B = B, s2 = s2, mu = mu, fy = fy, z = z,
             x.nds = x.nds, x.wts = x.wts,
             c.nds = c.nds, c.wts = c.wts,
             y.nds = y.nds, y.wts = y.wts)
  b <- get.b(B = B, s2 = s2, mu = mu, d.mu = d.mu,
             fy = fy, SF = SF, z = z,
             x.nds = x.nds, x.wts = x.wts,
             c.nds = c.nds, c.wts = c.wts,
             y.nds = y.nds, y.wts = y.wts)

  a <- t(chol2inv(diag(d) + M)) %*% t(b) / x.wts

  return(a)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' compute a function at nodes
#'
#' @inheritParams get.b
#'
#' @return approximated a values at nodes
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.a <- function(B, s2, mu, d.mu, fy, SF, zs,
                  x.nds, x.wts, c.nds, c.wts, y.wts, y.nds) {

  a <- lapply(X = 1:length(zs),
         FUN = function(i) get.a.z(
           B = B, s2 = s2, mu = mu, d.mu = d.mu, fy = fy, SF = SF, z = zs[i],
           x.nds = x.nds[,i], x.wts = x.wts[,i],
           c.nds = c.nds[,i], c.wts = c.wts[,i],
           y.nds = y.nds, y.wts = y.wts))

  return(a)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' interpolate a function between nodes
#'
#' @param a.vals a numeric matrix, values of a function at nodes
#' @param x.nds a numeric vector, nodes for quadrature rule for X
#' @param x.wts a numeric vector, weights for quadrature rule for X
#' @param x.new a numeric vector, x values at which to interpolate a
#' @param eta1 a function, density of X
#'
#' @return interpolated values of a
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
interp.a <- function(a.vals, x.nds, x.wts, x.new, z, eta1) {

  dim.a <- ncol(a.vals)
  dim.x <- length(x.new)
  a.new <- matrix(0, nrow = dim.x, ncol = dim.a)

  # loop through components of a() function
  for (i in 1:ncol(a.new)) {
    # interpolate a * x.wts
    a.new[,i] <- approx(x = x.nds,
                        y = a.vals[,i] * x.wts,
                        xout = x.new, rule = 2)$y
  }
  a.new <- a.new * sum(x.wts) / eta1(x.new, z)
  return(a.new)
}

