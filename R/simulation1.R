#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' simulation 1: assess various methods under misspecification
#'
#' @inheritParams get.q.beta
#' @inheritParams gen.data.beta
#'
#' @param B2 a number, the coefficient for `X` in the outcome model
#' @param s2 a number, the variance in the outcome model
#' @param mx a positive number, the number of nodes in quadrature grid for X
#' @param mc a positive number, the number of nodes in quadrature grid for C
#' @param my a positive number, the number of nodes in quadrature grid for Y
#' @param seed a positive integer, the seed value for random number generation
#'
#' @importFrom statmod gauss.quad
#' @importFrom fitdistrplus fitdistcens
#' @importFrom zipfR Ibeta
#' @import dplyr
#'
#' @export
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sim1 <- function(n, q, seed) {

  ## for troublehsooting
  #rm(list = ls()); library(devtools); load_all(); library(numDeriv)
  #n = 8000; q = 0.8; seed <- 1

  # define parameters -------------------------------------------------------

  # outcome model parameters
  B <- c(1, 10, 2)
  s2 <- 1

  # nuisance model parameters
  x.thetas <- 0.5 * c(-1, 1)
  x.gamma <- 1
  c.gamma <- 2

  # quadrature rule parameters
  mx <- 40
  mc <- 40
  my <- 5

  # spline estimation parameters
  m.knots <- 5
  deg <- 3

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

  ## X|Z correctly modeled as conditional beta
  st <- Sys.time()
  x.params.hat.correct <- t(vapply(
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
  et <- Sys.time()
  time.xz.param <- as.numeric(difftime(et, st, units = "secs"))

  eta1.correct <- function(x, z) {
    dbeta(x = x,
          shape1 = x.params.hat.correct[z + 1, "shape1"],
          shape2 = x.params.hat.correct[z + 1, "shape2"])
  }

  x.wts.correct <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mx),
    FUN = function(i) eta1.correct(x.nds[,i], zs[i]) /
      sum(eta1.correct(x.nds[,i], zs[i])))

  ## X|Z incorrectly modeled as marginal beta
  x.fit.wrong <- dat %>%
    mutate(left = W,
           right = ifelse(Delta == 1, W, NA)) %>%
    dplyr::select(left, right) %>%
    fitdistrplus::fitdistcens(distr = "beta")
  x.params.hat.wrong <- x.fit.wrong$estimate

  eta1.wrong <- function(x, z) {
    dbeta(x = x,
          shape1 = x.params.hat.wrong["shape1"],
          shape2 = x.params.hat.wrong["shape2"])
  }

  x.wts.wrong <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mx),
    FUN = function(i) eta1.wrong(x.nds[,i], zs[i]) /
      sum(eta1.wrong(x.nds[,i], zs[i])))

  ## nonparametric estimated distribution of X|Z (using B-splines)
  st <- Sys.time()
  spline.res.x <- fit.spline(dat = dat, m.knots = m.knots,
                             deg = deg, Boundary.knots = c(0, 1))
  eta1.nonpar <- spline.res.x$dens
  theta.hat.x <- spline.res.x$theta
  knots.x <- spline.res.x$knots
  wts.x <- spline.res.x$wts
  glogf <- spline.res.x$glogf
  et <- Sys.time()
  time.xz.nonpar <- as.numeric(difftime(et, st, units = "secs"))

  x.wts.nonpar <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mx),
    FUN = function(i) eta1.nonpar(x.nds[,i], zs[i]) /
      sum(eta1.nonpar(x.nds[,i], zs[i])))

  ## C|Z quadrature nodes
  c.nds <- vapply(
    X = zs,
    FUN.VALUE = numeric(mc),
    FUN = function(z) seq(1E-6, 1-1E-6, length = mc))

  ## C|Z correctly modeled as conditional beta
  st <- Sys.time()
  c.params.hat.correct <- t(vapply(
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
  et <- Sys.time()
  time.cz.param <- as.numeric(difftime(et, st, units = "secs"))

  eta2.correct <- function(c, z) {
    dbeta(x = c,
          shape1 = c.params.hat.correct[z + 1, "shape1"],
          shape2 = c.params.hat.correct[z + 1, "shape2"])
  }

  c.wts.correct <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mc),
    FUN = function(i) eta2.correct(c.nds[,i], zs[i]) /
      sum(eta2.correct(c.nds[,i], zs[i])))

  ## C|Z incorrectly modeled as marginal beta
  c.fit.wrong <- dat %>%
    mutate(left = W,
           right = ifelse(Delta == 0, W, NA)) %>%
    dplyr::select(left, right) %>%
    fitdistrplus::fitdistcens(distr = "beta")
  c.params.hat.wrong <- c.fit.wrong$estimate

  eta2.wrong <- function(c, z) {
    dbeta(x = c,
          shape1 = c.params.hat.wrong["shape1"],
          shape2 = c.params.hat.wrong["shape2"])
  }

  c.wts.wrong <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mc),
    FUN = function(i) eta2.wrong(c.nds[,i], zs[i]) /
      sum(eta2.wrong(c.nds[,i], zs[i])))

  ## nonparametric estimated distribution of C|Z (using B-splines)
  st <- Sys.time()
  spline.res.c <- fit.spline(dat = mutate(dat, Delta = 1 - Delta),
                             m.knots = m.knots, deg = deg,
                             Boundary.knots = c(0, 1))
  et <- Sys.time()
  time.cz.nonpar <- as.numeric(difftime(et, st, units = "secs"))

  eta2.nonpar <- spline.res.c$dens

  c.wts.nonpar <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mc),
    FUN = function(i) eta2.nonpar(c.nds[,i], zs[i]) /
      sum(eta2.nonpar(c.nds[,i], zs[i])))

  ## Y|X,Z quadrature nodes
  gq <- gauss.quad(n = my, kind = "hermite")
  y.nds <- gq$nodes
  y.wts <- gq$weights

  # estimate beta -----------------------------------------------------------

  ## complete case linear model to get starting value
  naive.lm <- lm(Y ~ W + Z, data = datcc)

  ## complete case
  st <- Sys.time()
  B.cc <- get.root(dat = dat, score = get.Scc,
                   start = c(naive.lm$coef, log(var(naive.lm$resid))))
  et <- Sys.time()
  time.est.cc <- as.numeric(difftime(et, st, units = "secs"))
  st <- Sys.time()
  V.cc <- var.est.sand(dat = datcc, theta = B.cc, args = list(),
                       n = sum(dat$Delta),
                       get.S = get.Scc, return.se = T)
  et <- Sys.time()
  time.var.cc <- as.numeric(difftime(et, st, units = "secs"))


  ## oracle
  st <- Sys.time()
  B.or <- get.root(dat = dat0, score = get.Scc, start = B.cc)
  et <- Sys.time()
  time.est.or <- as.numeric(difftime(et, st, units = "secs"))
  st <- Sys.time()
  V.or <- var.est.sand(dat = dat0, theta = B.or, args = list(),
                       get.S = get.Scc, return.se = T)
  et <- Sys.time()
  time.var.or <- as.numeric(difftime(et, st, units = "secs"))

  ## MLE (X|Z correct)
  mle.args <- list(mu = mu, d.mu = d.mu, SF = SF, fy = fy, x.nds = x.nds)
  st <- Sys.time()
  B.ml.1 <- get.root(dat = dat, score = get.Sml, start = B.cc,
                     args = append(mle.args, list("x.wts" = x.wts.correct)))
  et <- Sys.time()
  time.est.ml <- as.numeric(difftime(et, st, units = "secs"))
  st <- Sys.time()
  V.ml.1 <- tryCatch(
    expr =
      var.est.sand(
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
          S <- cbind(get.Sml(dat = dat, theta = head(theta, 4),
                             args = args, return.sums = F),
                     d.log.fxz(dat = dat, theta = alpha,
                               args = args, return.sums = F))

          if (return.sums) {
            return(colSums(S))
          } else {
            return(S)
          }
        },
        theta = c(B.ml.1, x.params.hat.correct),
        args = list(),
        return.se = T),
    error = rep(NA, length(V.cc)))
  et <- Sys.time()
  time.var.ml <- as.numeric(difftime(et, st, units = "secs"))

  ## MLE (X|Z incorrect)
  B.ml.0 <- get.root(dat = dat, score = get.Sml, start = B.cc,
                     args = append(mle.args, list("x.wts" = x.wts.wrong)))
  V.ml.0 <- tryCatch(
    expr =
    var.est.sand(
    dat = dat,
    get.S = function(dat, theta, args, return.sums = F) {

      alpha <- tail(theta, -4)

      # define estimated X|Z density
      eta1 <- function(x, z) {
        dbeta(x = x,
              shape1 = alpha[1],
              shape2 = alpha[2])
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
                 d.log.fx(dat = dat, theta = alpha,
                          args = args, return.sums = F))

      if (return.sums) {
        return(colSums(S))
      } else {
        return(S)
      }
    },
    theta = c(B.ml.0, x.params.hat.wrong),
    args = list(),
    return.se = T),
    error = rep(NA, length(V.cc)))

  ## MLE (X|Z nonparametric)
  B.ml.2 <- get.root(dat = dat, score = get.Sml, start = B.cc,
                     args = append(mle.args, list("x.wts" = x.wts.nonpar)))
  V.ml.2 <- var.est.sand(
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
          knots = knots.x,
          Boundary.knots = c(0, 1),
          degree = deg,
          intercept = F)

        # transform theta parameter to alpha
        alpha.hat <- apply(theta.nuis, 2, function(x) softmax(x, wts = wts.x))

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
                 glogf(dat = dat, theta = matrix(theta.nuis,
                                                 ncol = length(zs))))

      if (return.sums) {
        return(colSums(S))
      } else {
        return(S)
      }

    },
    theta = c(B.ml.2, theta.hat.x),
    args = list(),
    ridge.size = 1E-6,
    return.se = T
  )

  ## semiparametric (X|Z, C|Z correct)
  sp.args <- list(mu = mu, d.mu = d.mu, SF = SF, fy = fy,
                  x.nds = x.nds, c.nds = c.nds,
                  y.nds = y.nds, y.wts = y.wts)
  st <- Sys.time()
  B.sp.11 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.correct,
                                                  "c.wts" = c.wts.correct,
                                                  "eta1" = eta1.correct)))
  et <- Sys.time()
  time.est.sp <- as.numeric(difftime(et, st, units = "secs"))
  st <- Sys.time()
  V.sp.11 <- var.est.sand(dat = dat, theta = B.sp.11,
                          args = append(sp.args, list("x.wts" = x.wts.correct,
                                                      "c.wts" = c.wts.correct,
                                                      "eta1" = eta1.correct)),
                          get.S = get.Seff, return.se = T)
  et <- Sys.time()
  time.var.sp <- as.numeric(difftime(et, st, units = "secs"))

  ## semiparametric (X|Z correct, C|Z incorrect)
  B.sp.10 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.correct,
                                                  "c.wts" = c.wts.wrong,
                                                  "eta1" = eta1.correct)))
  V.sp.10 <- var.est.sand(dat = dat, theta = B.sp.10,
                          args = append(sp.args, list("x.wts" = x.wts.correct,
                                                      "c.wts" = c.wts.wrong,
                                                      "eta1" = eta1.correct)),
                          get.S = get.Seff, return.se = T)

  ## semiparametric (X|Z incorrect, C|Z correct)
  B.sp.01 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.wrong,
                                                  "c.wts" = c.wts.correct,
                                                  "eta1" = eta1.wrong)))
  V.sp.01 <- var.est.sand(dat = dat, theta = B.sp.01,
                          args = append(sp.args, list("x.wts" = x.wts.wrong,
                                                      "c.wts" = c.wts.correct,
                                                      "eta1" = eta1.wrong)),
                          get.S = get.Seff, return.se = T)

  ## semiparametric (X|Z, C|Z incorrect)
  B.sp.00 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.wrong,
                                                  "c.wts" = c.wts.wrong,
                                                  "eta1" = eta1.wrong)))
  V.sp.00 <- var.est.sand(dat = dat, theta = B.sp.00,
                          args = append(sp.args, list("x.wts" = x.wts.wrong,
                                                      "c.wts" = c.wts.wrong,
                                                      "eta1" = eta1.wrong)),
                          get.S = get.Seff, return.se = T)

  ## semiparametric (X|Z, C|Z nonparametric)
  B.sp.22 <- get.root(dat = dat, score = get.Seff, start = B.cc,
                      args = append(sp.args, list("x.wts" = x.wts.nonpar,
                                                  "c.wts" = c.wts.nonpar,
                                                  "eta1" = eta1.nonpar)))
  V.sp.22 <- var.est.sand(dat = dat, theta = B.sp.22,
                          args = append(sp.args, list("x.wts" = x.wts.nonpar,
                                                      "c.wts" = c.wts.nonpar,
                                                      "eta1" = eta1.nonpar)),
                          get.S = get.Seff, return.se = T)

  # return setup parameters, estimates, and standard errors (length 95)
  ret <- c(n = n, q = q, seed = seed,
           B.or = B.or, B.cc = B.cc,
           B.ml.0 = B.ml.0, B.ml.1 = B.ml.1, B.ml.2 = B.ml.2,
           B.sp.00 = B.sp.00, B.sp.01 = B.sp.01, B.sp.10 = B.sp.10,
           B.sp.11 = B.sp.11, B.sp.22 = B.sp.22,
           V.or = V.or, V.cc = V.cc,
           V.ml.0 = V.ml.0[1:4], V.ml.1 = V.ml.1[1:4], V.ml.2 = V.ml.2[1:4],
           V.sp.00 = V.sp.00, V.sp.01 = V.sp.01, V.sp.10 = V.sp.10,
           V.sp.11 = V.sp.11, V.sp.22 = V.sp.22,
           time.xz.param = time.xz.param, time.xz.nonpar = time.xz.nonpar,
           time.cz.param = time.cz.param, time.cz.nonpar = time.cz.nonpar,
           time.est.or = time.est.or, time.var.or = time.var.or,
           time.est.cc = time.est.cc, time.var.cc = time.var.cc,
           time.est.ml = time.est.ml, time.var.ml = time.var.ml,
           time.est.sp = time.est.sp, time.var.sp = time.var.sp)

  return(ret)
}


