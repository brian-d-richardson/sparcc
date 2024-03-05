#' Assess estimating equation
#'
#' @param ee a numeric matrix, estimating function values at each observation
#' @param digits a non-negative integer, number of digits to round results,
#' default is 2
#'
#' @return individual or summation estimating function values
#'
#' @export
assess.ee <- function(ee, digits = 2) {

  # means of estimating equation values (for each component)
  means <- colMeans(ee)

  # t statistics of H_0: mean = 0 (for each component)
  n <- nrow(ee)
  t.stats <- apply(ee, 2, function(x) mean(x) * sqrt(n) / sd(x))

  return(list(means = round(means, digits),
              t.stats = round(t.stats, digits)))
}


#' Assess data generation
#'
#' @inheritParams gen.data
#' @inheritParams get.q
#'
#' @param x.mean a positive number, the mean of X
#'
#' @importFrom ggplot2 ggplot geom_histogram labs ggtitle
#'
#' @return a plot of generated data
#'
#' @export
assess.dat <- function(n, q, B, s2, x.mean, x.shape, c.shape,
                       specify.x.gamma, specify.c.gamma) {

  x.rate <- x.shape / x.mean # rate parameter for gamma distribution of X
  c.rate <- get.c.rate(      # rate parameter for gamma distribution of C
    q = q,
    x.rate = x.rate,
    x.shape = x.shape,
    c.shape = c.shape)

  dat.list <- gen.data(n = n, q = q, B = B, s2 = s2,
                       x.rate = x.rate, x.shape = x.shape,
                       c.rate = c.rate, c.shape = c.shape)
  datf <- dat.list$datf          # full data
  dat0 <- dat.list$dat0          # oracle data
  dat <- dat.list$dat            # observed data
  datcc <- dat.list$datcc        # complete case data

  # compute empirical values
  q.hat <- mean(datf$X > datf$C)                   # censoring proportion
  x.bar <- mean(datf$X)                            # sample mean of X

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

  x.grid <- seq(0, max(c(datf$X, datf$C)), 0.01)   # grid to plot X density

  # plot data
  ggplot() +
    geom_histogram(data = datf,
                   aes(x = X,
                       y = after_stat(density)),
                   fill = "blue", alpha = 0.5, bins = n/50) +
    geom_histogram(data = datf,
                   aes(x = C,
                       y = after_stat(density)),
                   fill = "red", alpha = 0.5, bins = n/50) +
    geom_line(aes(x = x.grid,
                  y = eta1(x.grid)),
              color = "blue", linewidth = 1) +
    geom_line(aes(x = x.grid,
                  y = eta2(x.grid)),
              color = "red", linewidth = 1) +
    labs(x = "X (blue) or C (red)",
         y = "Count") +
    ggtitle("Estimated vs Observed Distributions of X and C",
            subtitle = paste0("q.hat = ", round(q.hat, 2), "; ",
                              "x.bar = ", round(x.bar, 2)))
}


#' Assess estimation of eta
#'
#' @inheritParams gen.data
#' @inheritParams get.q
#'
#' @param x.mean a positive number, the mean of X
#' @param n.rep a positive integer, the number of simulated replicates
#'
#' @return empirical distribution of estimated parameters for eta1 and eta2
#'
#' @export
assess.eta.est <- function(n.rep, n, q, B, s2, x.mean, x.shape, c.shape,
                           specify.x.gamma, specify.c.gamma) {

  # rate parameters for X and C
  x.rate <- x.shape / x.mean
  c.rate <- get.c.rate(
    q = q,
    x.rate = x.rate,
    x.shape = x.shape,
    c.shape = c.shape)

  param.hat <- t(pbapply::pbvapply(
    X = 1:n.rep,
    FUN.VAL = numeric(4),
    FUN = function(x) {

        # generate data
        dat <- gen.data(n = n, q = q, B = B, s2 = s2,
                        x.rate = x.rate, x.shape = x.shape,
                        c.rate = c.rate, c.shape = c.shape)$dat

        # estimate distribution of X
        if (specify.x.gamma) {
          x.param.hat <- gammaMLE(yi = dat$W,                # gamma params
                                  si = dat$Delta,
                                  scale = F)$estimate
        } else {
          x.param.hat <- c(1, mean(dat$Delta) / mean(dat$W)) # exp params
        }

        # estimate distribution of C
        if (specify.c.gamma) {
          c.param.hat <- gammaMLE(yi = dat$W,                    # gamma params
                                  si = 1 - dat$Delta,
                                  scale = F)$estimate
        } else {
          c.param.hat <- c(1, mean(1 - dat$Delta) / mean(dat$W)) # exp params
        }
        ret <- c(x.param.hat, c.param.hat)
        names(ret) <- c("x.shape", "x.rate", "c.shape", "c.rate")
        return(ret)
      })) |>
    as.data.frame()

  plot.dat <- data.frame(
    var = rep(c("X", "C"), each = 2 * n.rep),
    param = rep(c("shape", "rate", "shape", "rate"), each = n.rep),
    hat = c(param.hat$x.shape, param.hat$x.rate,
            param.hat$c.shape, param.hat$c.rate),
    true = rep(c(x.shape, x.rate, c.shape, c.rate), each = n.rep))

  plot <- ggplot(
    data = plot.dat,
    aes(y = hat,
        fill = param)) +
    geom_boxplot() +
    geom_hline(aes(yintercept = true),
               linetype = "dashed",
               color = "orange") +
    scale_fill_manual(values = c("red", "blue")) +
    facet_wrap(var ~ param, scales = "free_y") +
    labs(y = "Estimated Parameter") +
    ggtitle("Empirical Distribution of Estimated Nuisance Parameters",
            subtitle = paste0("q = ", q, "; ",
                              "n = ", n, "; ",
                              "n.rep = ", n.rep, "\n",
                              "x.shape = ", x.shape, "; ",
                              "c.shape = ", c.shape)) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none")

  return(plot)
}
