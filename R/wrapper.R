#' Fit SPARCC Estimator
#'
#' @description
#' To fit the SPARCC estimator for linear model with a conditionally normal
#' outcome `Y`, a censored covariate `X`, and a discrete-valued uncensored
#' covariate `Z`.
#'
#'
#' @param data a data frame with the following columns
#' \itemize{
#' \item{`Y`: outcome}
#' \item{`W`: minimum of the censored covariate `X` and the censoring time `C`}
#' \item{`Delta`: an indicator if `X` is observed (`Delta = 1`) or censored
#' (`Delta = 0`)}
#' \item{`Z`: a discrete-valued uncensored covariate}
#' }
#'
#' @param xz.interaction a logical indicator for whether an `X*Z` interaction is
#' to be included in the outcome model (`TRUE`) or not (`FALSE`), default is
#' `TRUE`
#'
#' @param nuisance.models model type for nuisance distribution models
#' ("parametric" or "nonparametric)
#'
#' @param distr.x a character string "name" naming a distribution for `X|Z`,
#' for which the corresponding density function dname must be defined; used only
#' if `nuisance.models` is "parametric"
#'
#' @param distr.c a character string "name" naming a distribution for `C|Z`,
#' for which the corresponding density function dname must be defined; used only
#' if `nuisance.models` is "parametric"
#'
#' @param m.knots a positive integer, the number of knots in the nonparametric
#' B-spline estimator for the densities of of `X|Z` and `C|Z`, default is 5;
#' used only if `nuisance.models` is "nonparametric"
#'
#' @param deg a positive integer, the degree of polynomials in the nonparametric
#' B-spline estimator for the densities of of `X|Z` and `C|Z`, default is 3;
#' used only if `nuisance.models` is "nonparametric"
#'
#' @param mx a positive integer, the number of quadrature nodes used for `X`,
#' default is 40
#'
#' @param mc a positive integer, the number of quadrature nodes used for `C`,
#' default is 40
#'
#' @param my a positive integer, the number of Gauss-Hermite quadrature nodes
#' used for `Y`, default is 5
#'
#' @param range.x a numeric vector of length 2, the range of the support of `X`,
#' default is c(1E-6, 1-1E-6)
#'
#' @param range.c a numeric vector of length 2, the range of the support of `C`,
#' default is c(1E-6, 1-1E-6)
#'
#' @return a list of arguments including
#' \itemize{
#' \item{`x.model`: a list with
#' `x.model.type` (parametric or nonparametric),
#' `distr.x`, `m.knots`, `deg`,
#' `x.params.hat` (the estimated parameters for the distribution of `X|Z` at
#' each level of `Z`),
#' `x.nds` (the quadrature nodes used for `X`),
#' `x.wts` (the quadrature weights used for `X`),
#' `eta1` (the estimated density function of `X|Z` at specific nodes)
#' }
#' \item{`c.model`: a list with
#' `c.model.type` (parametric or nonparametric),
#' `distr.c`, `m.knots`, `deg`,
#' `c.params.hat` (the estimated parameters for the distribution of `C|Z` at
#' each level of `Z`),
#' `c.nds` (the quadrature nodes used for `C`),
#' `c.wts` (the quadrature weights used for `C`),
#' `eta2` (the estimated density function of `C|Z` at specific nodes)
#' }
#' \item{`outcome.model`: a list with
#' `outcome.fmla` (the outcome model formula),
#' `coef` (the outcome model coefficient estimates),
#' `cov` (the estimated covariance matrix of the coefficient estimator)}
#' }
#'
#' @export
sparcc <- function(

    data,

    # for outcome model speification
    xz.interaction = T,

    # for nuisance model estimation
    nuisance.models = "parametric",

    # for parametric nuisance model estimation
    distr.x = "beta",
    distr.c = "beta",

    # for nonparametric nuisance model estimation
    m.knots = 5,
    deg = 3,

    # for creating quadrature rules:
    mx = 40,
    my = 40,
    mc = 5,
    range.x = c(1E-6, 1-1E-6),
    range.c = c(1E-6, 1-1E-6)

    ) {

  # for troubleshooting
  #data = dat
  #nuisance.models = "parametric"; distr.x = "beta"; distr.c = "beta"
  #mx = 40; mc = 40;  my = 5
  #range.x = c(1E-6, 1-1E-6); range.c = c(1E-6, 1-1E-6)

  # extract values ----------------------------------------------------------

  ## unique z values
  zs <- sort(unique(dat$Z))

  ## outcome model
  outcome.fmla <- if (xz.interaction) {
    Y ~ X * Z
  } else {
    Y ~ X + Z
  }

  ## X|Z quadrature nodes
  x.nds <- vapply(
    X = zs,
    FUN.VALUE = numeric(mx),
    FUN = function(z) seq(range.x[1], range.x[2], length = mx))

  ## C|Z quadrature nodes
  c.nds <- vapply(
    X = zs,
    FUN.VALUE = numeric(mc),
    FUN = function(z) seq(range.c[1], range.c[2], length = mc))

  ## Y|X,Z gauss-hermite quadrature
  gq <- gauss.quad(n = my, kind = "hermite")
  y.nds <- gq$nodes
  y.wts <- gq$weights

  # fit nuisance models -----------------------------------------------------

  message(paste0("STEP 1: fit ", nuisance.models, " nuisance models"))
  st1 <- Sys.time()

  if (nuisance.models == "parametric") {

    ## density function for distributions of X|Z and C|Z
    dens.x <- get(paste0("d", distr.x))
    dens.c <- get(paste0("d", distr.c))

    ## fit parametric model for X|Z
    x.params.hat <- t(vapply(
      X = 0:1,
      FUN.VALUE = numeric(2),
      FUN = function(z) {
        est <- data %>%
          filter(Z == z) %>%
          mutate(left = W,
                 right = ifelse(Delta == 1, W, NA)) %>%
          dplyr::select(left, right) %>%
          { suppressWarnings(fitdistrplus::fitdistcens(., distr = distr.x)) }
        return(est$estimate)
      }))
    eta1.hat <- function(x, z) {
      do.call(dens.x, c(list(x = x), x.params.hat[z + 1, ]))
    }

    ## fit parametric model for C|Z
    c.params.hat <- t(vapply(
      X = 0:1,
      FUN.VALUE = numeric(2),
      FUN = function(z) {
        est <- dat %>%
          filter(Z == z) %>%
          mutate(left = W,
                 right = ifelse(Delta == 0, W, NA)) %>%
          dplyr::select(left, right) %>%
          { suppressWarnings(fitdistrplus::fitdistcens(., distr = distr.c)) }
        return(est$estimate)
      }))
    eta2.hat <- function(c, z) {
      do.call(dens.c, c(list(x = c), c.params.hat[z + 1, ]))
    }

    ## nullify nonparametric model values
    m.knots <- NULL
    deg <- NULL

  } else if (nuisance.models == "nonparametric") {

    ## boundary knots for X|Z and C|Z
    bdry.x <- range.x + 1E-6 * c(-1, 1)
    bdry.c <- range.c + 1E-6 * c(-1, 1)

    ## fit nonparametric model for X|Z
    spline.res.x <- fit.spline(dat = data, m.knots = m.knots,
                               deg = deg, Boundary.knots = bdry.x)
    eta1.hat <- spline.res.x$dens
    theta.hat.x <- spline.res.x$theta
    knots.x <- spline.res.x$knots
    wts.x <- spline.res.x$wts
    glogf.x <- spline.res.x$glogf

    ## fit nonparametric model for C|Z
    spline.res.c <- fit.spline(dat = mutate(data, Delta = 1 - Delta),
                               m.knots = m.knots, deg = deg,
                               Boundary.knots = bdry.c)
    eta2.hat <- spline.res.c$dens

    ## nullify parametric model values
    distr.x = NULL
    distr.c = NULL
    x.params.hat = NULL
    c.params.hat = NULL
  }

  ## create X|Z quadrature weights
  x.wts <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mx),
    FUN = function(i) eta1.hat(x.nds[,i], zs[i]) /
      sum(eta1.hat(x.nds[,i], zs[i])))

  ## create C|Z quadrature weights
  c.wts <- vapply(
    X = 1:length(zs),
    FUN.VALUE = numeric(mc),
    FUN = function(i) eta2.hat(c.nds[,i], zs[i]) /
      sum(eta2.hat(c.nds[,i], zs[i])))

  et1 <- Sys.time()
  t1 <- round(as.numeric(difftime(et1, st1, units = "secs")), 2)
  message(paste0("STEP 1 complete (", t1, " seconds)"))

  # fit outcome model -------------------------------------------------------

  message(paste0("STEP 2: obtain SPARCC estimator"))
  st2 <- Sys.time()

  ## complete case linear model to get starting value
  naive.fmla <- if (xz.interaction) {
    Y ~ W * Z
  } else {
    Y ~ W + Z
  }
  cc.lm <- lm(naive.fmla,
              data = filter(data, Delta == 1))
  start <- c(cc.lm$coef, log(var(cc.lm$resid)))

  sp.args <- list(
    mu = mu, d.mu = d.mu, SF = SF, fy = fy,
    eta1 = eta1.hat,
    xz.interaction = xz.interaction,
    x.nds = x.nds, x.wts = x.wts,
    c.nds = c.nds, c.wts = c.wts,
    y.nds = y.nds, y.wts = y.wts)

  sparcc.root <- get.root(
    dat = data, score = get.Seff,
    args = sp.args, start = start)

  et2 <- Sys.time()
  t2 <- round(as.numeric(difftime(et2, st2, units = "secs")), 2)
  message(paste0("STEP 2 complete (", t2, " seconds)"))

  # variance estimation -----------------------------------------------------

  message(paste0("STEP 3: obtain SPARCC variance estimator"))
  st3 <- Sys.time()

  sparcc.var <- var.est.sand(
    dat = data,
    theta = sparcc.root,
    args = sp.args,
    get.S = get.Seff,
    return.se = F)

  et3 <- Sys.time()
  t3 <- round(as.numeric(difftime(et3, st3, units = "secs")), 2)
  message(paste0("STEP 3 complete (", t3, " seconds)"))


  # format quadrature rules for output --------------------------------------

  ## X|Z quadrature
  eta1 <- cbind(
    x.nds,
    vapply(
      X = 1:length(zs),
      FUN.VALUE = numeric(mx),
      FUN = function(i) eta1.hat(x.nds[,i], zs[i]))) %>%
    `colnames<-`(do.call(paste0, expand.grid(zs, c("x.nds", "eta1")))) %>%
    as.data.frame() %>%
    mutate(row = row_number()) %>%
    tidyr::pivot_longer(
      cols = -row,
      names_to = c("Z", "var"),
      names_pattern = "([01])(.*)",
      values_to = "value") %>%
    tidyr::pivot_wider(
      id_cols = c(row, Z),
      names_from = var,
      values_from = value) %>%
    mutate(Z = as.integer(Z)) %>%
    select(Z, x.nds, eta1)

  ## C|Z quadrature
  eta2 <- cbind(
    c.nds,
    vapply(
      X = 1:length(zs),
      FUN.VALUE = numeric(mc),
      FUN = function(i) eta2.hat(c.nds[,i], zs[i]))) %>%
    `colnames<-`(do.call(paste0, expand.grid(zs, c("c.nds", "eta2")))) %>%
    as.data.frame() %>%
    mutate(row = row_number()) %>%
    tidyr::pivot_longer(
      cols = -row,
      names_to = c("Z", "var"),
      names_pattern = "([01])(.*)",
      values_to = "value") %>%
    tidyr::pivot_wider(
      id_cols = c(row, Z),
      names_from = var,
      values_from = value) %>%
    mutate(Z = as.integer(Z)) %>%
    select(Z, c.nds, eta2)

  # format output -----------------------------------------------------------

  ## X|Z model
  x.model <- list(
    x.model.type = nuisance.models,
    distr.x = distr.x,
    m.knots = m.knots,
    deg = deg,
    x.params.hat = x.params.hat,
    eta1 = eta1
  )

  ## C|Z model
  c.model <- list(
    c.model.type = nuisance.models,
    distr.c = distr.c,
    m.knots = m.knots,
    deg = deg,
    c.params.hat = c.params.hat,
    eta2 = eta2
  )

  ## outcome model
  outcome.model <- list(
    outcome.fmla = outcome.fmla,
    coef = sparcc.root,
    cov = sparcc.var
  )

  return(list(
    x.model = x.model,
    c.model = c.model,
    outcome.model = outcome.model
  ))

}
