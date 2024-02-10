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
assess.dat <- function(n, q, B, s2, x.mean, x.shape, c.shape) {

  x.rate <- x.shape / x.mean

  c.rate <- get.c.rate(      # rate parameter for gamma distribution of C
    q = q,
    x.rate = x.rate,
    x.shape = x.shape,
    c.shape = c.shape)

  # generate data
  dat <- gen.data(n = n, q = q, B = B, s2 = s2,
                  x.rate = x.rate, x.shape = x.shape,
                  c.rate = c.rate, c.shape = c.shape)
  datf <- dat$datf # full data (Y, X, C)
  dat <- dat$dat   # observed data (Y, W, Delta)

  # compute empirical values
  q.hat <- mean(datf$X > datf$C)                   # censoring proportion
  x.bar <- mean(datf$X)                            # sample mean of X
  x.rate.hat <- mean(dat$Delta) / mean(dat$W)      # mle for exponential X rate
  c.rate.hat <- mean(1 - dat$Delta) / mean(dat$W)  # mle for exponential C rate
  x.grid <- seq(0, max(c(datf$X, datf$C)), 0.01)

  # plot data
  ggplot() +
    geom_histogram(data = datf,
                   aes(x = X,
                       y = after_stat(density)),
                   fill = "blue", alpha = 0.5, bins = 200) +
    geom_histogram(data = datf,
                   aes(x = C,
                       y = after_stat(density)),
                   fill = "red", alpha = 0.5, bins = 200) +
    geom_line(aes(x = x.grid,
                  y = dexp(x = x.grid, rate = x.rate.hat)),
              color = "blue", linewidth = 1) +
    geom_line(aes(x = x.grid,
                  y = dexp(x = x.grid, rate = c.rate.hat)),
              color = "red", linewidth = 1) +
    labs(x = "X (blue) or C (red)",
         y = "Count") +
    ggtitle("Estimated vs Observed Distributions of X and C",
            subtitle = paste0("x.shape = ", x.shape, "; ",
                              "c.shape = ", c.shape, "\n",
                              "q.hat = ", round(q.hat, 2), "; ",
                              "x.bar = ", round(x.bar, 2), "\n",
                              "x.rate.hat = ", round(x.rate.hat, 2), "; ",
                              "c.rate.hat = ", round(c.rate.hat, 2)))
}
