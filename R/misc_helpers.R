#' Assess estimating equation
#'
#' @param ee ...
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
#' @param
#'
#' @return
#'
#' @importFrom ggplot2 ggplot geom_histogram labs ggtitle
#'
#' @export
assess.dat <- function(n, q, B, x.mean, x.shape, c.shape) {

  x.rate <- x.shape / x.mean

  c.rate <- get.c.rate(      # rate parameter for gamma distribution of C
    q = q,
    x.rate = x.rate,
    x.shape = x.shape,
    c.shape = c.shape)

  # generate full data (Y, X, C)
  datf <- gen.data(n = n, q = q, B = B,
                   x.rate = x.rate, x.shape = x.shape,
                   c.rate = c.rate, c.shape = c.shape)$datf

  # compute empirical values
  q.hat <- mean(datf$X > datf$C)  # censoring proportion
  x.bar <- round(mean(datf$X), 2) # sample mean of X

  # plot data
  ggplot(data = datf) +
    geom_histogram(aes(x = X), fill = "blue", alpha = 0.5, bins = 30) +
    geom_histogram(aes(x = C), fill = "red", alpha = 0.5, bins = 30) +
    labs(x = "X (blue) or C (red)",
         y = "Count") +
    ggtitle("Distributions of X and C",
            subtitle = paste0("x.shape = ", x.shape, "; ",
                              "c.shape = ", c.shape, "; ",
                              "q.hat = ", q.hat, "; ",
                              "x.bar = ", x.bar))
}
