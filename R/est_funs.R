#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' complete case score function
#'
#' @param dat a data frame including columns Y, W, Delta, Z
#' @param B a numeric vector, parameters in the outcome model
#' @param s2 a positive number, variance in the outcome model
#' @param args list of additional arguments
#' @param return.sums logical indicator for returning sum of scores as opposed
#' to individual scores, default is TRUE
#'
#' @return complete case score function values
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.Scc <- function(dat, B, ls2, args, return.sums = T) {

  # unpack arguments
  list2env(args, envir = environment())

  # full score for uncensored observations
  Scc <- SF(
    y = dat$Y[dat$Delta == 1],
    x = dat$W[dat$Delta == 1],
    z = dat$Z[dat$Delta == 1],
    B = B, ls2 = ls2)

  if (return.sums) {
    return(colSums(Scc))
  } else {
    return(Scc)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' MLE score function
#'
#' @inheritParams get.Scc
#'
#' @return maximum likelihood score function values
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.Sml <- function(dat, B, ls2, args, return.sums = T) {

  # unpack arguments
  list2env(args, envir = environment())

  # unique levels of z
  zs <- sort(unique(dat$Z))

  # initiate score vector
  Sml <- matrix(nrow = nrow(dat), ncol = length(B) + 1)

  # full score for uncensored observations
  Sml[dat$Delta == 1, ] <- SF(
    y = dat$Y[dat$Delta == 1],
    x = dat$W[dat$Delta == 1],
    z = dat$Z[dat$Delta == 1],
    B = B, ls2 = ls2)

  # expected score for censored observations
  for (i in which(dat$Delta == 0)) {

    # index for Z
    zi <- which(dat$Z[i] == zs)

    # x node indices greater than observed W = C
    xi <- x.nds[,zi] > dat$W[i]

    # (proportional to) joint density of Y, X
    fyx <- fy(y = dat$Y[i], x = x.nds[xi, zi],
              z = dat$Z[i], B = B, s2 = exp(ls2)) *
      x.wts[xi, zi]

    # conditional expectation of SF
    Sml[i,] <- colSums(SF(y = dat$Y[i], x = x.nds[xi, zi],
                          z = dat$Z[i], B = B, ls2 = ls2) * fyx) /
      sum(fyx)
  }

  if (return.sums) {
    return(colSums(Sml))
  } else {
    return(Sml)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' semiparametric efficient score function
#'
#' @inheritParams get.Scc
#'
#' @return semiparametric efficient score function values
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.Seff <- function(dat, B, ls2, args, return.sums = T) {

  # unpack arguments
  list2env(args, envir = environment())

  # unique levels of z
  zs <- sort(unique(dat$Z))

  # initiate score vector
  Seff <- matrix(nrow = nrow(dat), ncol = length(B) + 1)

  # solve integral equation for a() function at x nodes
  a.vals <- get.a(B = B, s2 = exp(ls2), mu = mu, d.mu = d.mu,
                  SF = SF, fy = fy, zs = zs,
                  x.nds = x.nds, x.wts = x.wts,
                  c.nds = c.nds, c.wts = c.wts,
                  y.nds = y.nds, y.wts = y.wts)

  # full score minus interpolated a() for uncensored observations
  for (z in zs) {
    zi <- which(z == zs)
    Seff[dat$Delta == 1 & dat$Z == z] <- SF(
      y = dat$Y[dat$Delta == 1 & dat$Z == z],
      x = dat$W[dat$Delta == 1 & dat$Z == z],
      z = z, B = B, ls2 = ls2) -
      interp.a(a.vals = a.vals[[zi]], x.nds = x.nds[,zi], x.wts = x.wts[,zi],
               x.new = dat$W[dat$Delta == 1 & dat$Z == z],
               z = z, eta1 = eta1)
  }

  # expected score minus a() for censored observations
  for (i in which(dat$Delta == 0)) {

    # index for Z
    zi <- which(dat$Z[i] == zs)

    # x node indices greater than observed W = C
    xi <- x.nds[,zi] > dat$W[i]

    # (proportional to) joint density of Y, X
    fyx <- fy(y = dat$Y[i], x = x.nds[xi, zi],
              z = dat$Z[i], B = B, s2 = exp(ls2)) *
      x.wts[xi, zi]

    # conditional expectation of SF
    Seff[i,] <- colSums(SF(y = dat$Y[i], x = x.nds[xi, zi],
                           z = dat$Z[i], B = B, ls2 = ls2) *
                          a.vals[[zi]][xi,] * fyx) /
      sum(fyx)
  }

  if (return.sums) {
    return(colSums(Seff))
  } else {
    return(Seff)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' find root of estimating function
#'
#' @inheritParams get.Scc
#' @param score a function of theta = c(beta, log(s2)), the estimating function
#' to be solved
#' @param start a numeric vector, starting value for root search
#'
#' @importFrom rootSolve multiroot
#'
#' @return root of estimating function
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.root <- function(dat, score, start, args = list()) {

  est <- #tryCatch(
    #expr =
    rootSolve::multiroot(
      f = function(theta) score(dat = dat, args = args,
                                B = head(theta, -1), ls2 = tail(theta, 1)),
      start = start)$root#,
    #warning = function(w) rep(NA, length(start)),
    #error = function(e) rep(NA, length(start)))
  names(est) <- paste0("B", 1:length(start))
  return(est)
}

