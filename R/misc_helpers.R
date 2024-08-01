#' Assess estimating equation
#'
#' @param ee a numeric matrix, estimating function values at each observation
#' @param digits a non-negative integer, number of digits to round results,
#' default is 2
#'
#' @return A list with:
#' #' \itemize{
#' \item{`means`: the mean estimatin equation values}
#' \item{`t.stats`: t-statistics testing the null of mean zero}
#' }
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


#' Assess data generation with X|Z and C|Z having gamma distributions
#'
#' @inheritParams gen.data.gamma
#' @inheritParams get.q.gamma
#'
#' @param x.mean a positive number, the mean of X
#'
#' @importFrom ggplot2 ggplot geom_histogram labs ggtitle
#'
#' @return a plot of generated data
#'
#' @export
assess.dat.gamma <- function(n, q, B, s2, x.means, x.shape, c.shape,
                             specify.x.gamma, specify.c.gamma) {

  dat.list <- gen.data.gamma(n = n, q = q, B = B, s2 = s2,
                             x.means = x.means,
                             x.shape = x.shape, c.shape = c.shape)

  datf <- dat.list$datf          # full data
  dat0 <- dat.list$dat0          # oracle data
  dat <- dat.list$dat            # observed data
  datcc <- dat.list$datcc        # complete case data

  # compute empirical values by Z level
  q.hats <- datf %>% group_by(Z) %>% summarise(qhat = mean(X > C)) # cens prop
  x.bars <- datf %>% group_by(Z) %>% summarise(xbar = mean(X))     # mean X

  # estimate distribution of X|Z
  if (specify.x.gamma) {

    # estimate gamma parameters at each level of Z
    x.params.hat <- t(vapply(
      X = 0:1,
      FUN.VALUE = numeric(2),
      FUN = function(z) {
        z.ind <- dat$Z == z
        gammaMLE(yi = dat$W[z.ind],
                 si = dat$Delta[z.ind],
                 scale = F)$estimate
      }))

    # define estimated X|Z density
    eta1 <- function(x, z) {
      dgamma(x = x,
             shape = x.params.hat[z + 1, "shape"],
             rate = x.params.hat[z + 1, "rate"])
    }

  } else {

    # estimate exponential parameter at each level of Z
    x.rates.hat <- vapply(
      X = 0:1,
      FUN.VALUE = 0,
      FUN = function(z) {
        z.ind <- dat$Z == z
        mean(dat$Delta[z.ind]) / mean(dat$W[z.ind])
      })

    # define estimated X|Z density
    eta1 <- function(x, z) {
      dexp(x, rate = x.rates.hat[z + 1])
    }
  }

  # estimate distribution of C|Z
  if (specify.c.gamma) {

    # estimate gamma parameters at each level of Z
    c.params.hat <- t(vapply(
      X = 0:1,
      FUN.VALUE = numeric(2),
      FUN = function(z) {
        z.ind <- dat$Z == z
        gammaMLE(yi = dat$W[z.ind],
                 si = 1 - dat$Delta[z.ind],
                 scale = F)$estimate
      }))

    # define estimated C|Z density
    eta2 <- function(c, z) {
      dgamma(x = c,
             shape = c.params.hat[z + 1, "shape"],
             rate = c.params.hat[z + 1, "rate"])
    }

  } else {

    # estimate exponential parameter at each level of Z
    c.rates.hat <- vapply(
      X = 0:1,
      FUN.VALUE = 0,
      FUN = function(z) {
        z.ind <- dat$Z == z
        mean(1 - dat$Delta[z.ind]) / mean(dat$W[z.ind])
      })

    # define estimated C|Z density
    eta2 <- function(c, z) {
      dexp(x = c, rate = c.rates.hat[z + 1])
    }
  }

  # grid to plot X density
  x.grid <- seq(0, max(c(datf$X, datf$C)), 0.01)

  # plot data
  datf %>%
    mutate(e1 = eta1(x = X, z = Z),
           e2 = eta2(c = C, z = Z)) %>%
    ggplot() +
    geom_histogram(aes(x = X,
                       y = after_stat(density)),
                   fill = "blue", alpha = 0.5, bins = n/50) +
    geom_histogram(aes(x = C,
                       y = after_stat(density)),
                   fill = "red", alpha = 0.5, bins = n/50) +
    geom_line(aes(x = X,
                  y = e1),
              color = "blue", linewidth = 1) +
    geom_line(aes(x = C,
                  y = e2),
              color = "red", linewidth = 1) +
    facet_wrap(~ Z,
               scales = "free", labeller = label_both) +
    labs(x = "X (blue) or C (red)",
         y = "Density") +
    ggtitle("Estimated vs Observed Distributions of X|Z and C|Z",
            subtitle = paste0("q.hats = ", paste(round(q.hats$qhat, 2), collapse = ", "), "; ",
                              "x.bars = ", paste(round(x.bars$xbar, 2), collapse = ", ")))

}

#' Assess estimation of eta when X and C have gamma distributions
#'
#' @inheritParams gen.data.gamma
#' @inheritParams get.q.gamma
#'
#' @param x.mean a positive number, the mean of X
#' @param n.rep a positive integer, the number of simulated replicates
#' @param specify.x.gamma logical, an indicator of whether X should be correctly
#' specified as gamma as opposed to incorrectly specified as exponential
#' @param specify.c.gamma logical, an indicator of whether C should be correctly
#' specified as gamma as opposed to incorrectly specified as exponential
#'
#' @return plot of empirical distribution of estimated parameters for eta1 and
#' eta2
#'
#' @importFrom tidyr separate_wider_delim
#' @importFrom ggh4x facet_nested
#'
#' @export
assess.eta.est.gamma <- function(n.rep, n, q, B, s2, x.means, x.shape, c.shape,
                                 specify.x.gamma, specify.c.gamma) {

  x.rates <- x.shape / x.means  # rate parameters for gamma distribution of X|Z
  c.rates <- vapply(            # rate parameters for gamma distribution of C|Z
    X = x.rates,
    FUN.VALUE = 0,
    FUN = function(xr)
      get.c.rate(q = q, x.rate = xr,
                 x.shape = x.shape,
                 c.shape = c.shape))

  param.hat <- t(pbapply::pbvapply(
    X = 1:n.rep,
    FUN.VAL = numeric(8),
    FUN = function(x) {

        # generate data
        dat <- gen.data.gamma(n = n, q = q, B = B, s2 = s2,
                              x.means = x.means,
                              x.shape = x.shape, c.shape = c.shape)$dat

        # estimate distribution of X
        if (specify.x.gamma) {

          # estimate gamma parameters at each level of Z
          x.params.hat <- t(vapply(
            X = 0:1,
            FUN.VALUE = numeric(2),
            FUN = function(z) {
              z.ind <- dat$Z == z
              gammaMLE(yi = dat$W[z.ind],
                       si = dat$Delta[z.ind],
                       scale = F)$estimate
            }))

        } else {

          # estimate exponential parameter at each level of Z
          x.params.hat <- vapply(
            X = 0:1,
            FUN.VALUE = numeric(2),
            FUN = function(z) {
              z.ind <- dat$Z == z
              c(1, mean(dat$Delta[z.ind]) / mean(dat$W[z.ind]))
            })
        }

        # estimate distribution of C
        if (specify.c.gamma) {

          # estimate gamma parameters at each level of Z
          c.params.hat <- t(vapply(
            X = 0:1,
            FUN.VALUE = numeric(2),
            FUN = function(z) {
              z.ind <- dat$Z == z
              gammaMLE(yi = dat$W[z.ind],
                       si = 1 - dat$Delta[z.ind],
                       scale = F)$estimate
            }))

        } else {

          # estimate exponential parameter at each level of Z
          c.params.hat <- vapply(
            X = 0:1,
            FUN.VALUE = 0,
            FUN = function(z) {
              z.ind <- dat$Z == z
              mean(1 - dat$Delta[z.ind]) / mean(dat$W[z.ind])
            })

        }
        ret <- c(x.params.hat, c.params.hat)
        #colnames(ret) <- c("z", "x.shape", "x.rate", "c.shape", "c.rate")
        return(ret)
      }))

  plot.dat <- cbind(rep(c(0, 1), each = n.rep),
                    rbind(param.hat[,c(1,3,5,7)],
                          param.hat[,c(2,4,6,8)])) %>%
    as.data.frame() %>%
    `colnames<-`(c("z", "x.shape", "x.rate", "c.shape", "c.rate")) %>%
    tidyr::pivot_longer(cols = !z, values_to = "est") %>%
    tidyr::separate_wider_delim(cols = name, delim = ".",
                                names = c("var", "param")) %>%
    mutate(truth = ifelse(var == "x",
                          ifelse(param == "rate", x.rates[z + 1], x.shape),
                          ifelse(param == "rate", c.rates[z + 1], c.shape)))

  plot <- ggplot(
    data = plot.dat,
    aes(y = est,
        fill = param)) +
    geom_boxplot() +
    geom_hline(aes(yintercept = truth),
               linetype = "dashed",
               color = "orange") +
    scale_fill_manual(values = c("red", "blue")) +
    ggh4x::facet_nested(param ~ var + z,
                        scales = "free",
                        labeller = label_both) +
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

#' Assess data generation with X|Z and C|Z having beta distributions
#'
#' @inheritParams gen.data.beta
#' @inheritParams get.q.beta
#'
#' @param x.mean a positive number, the mean of X
#'
#' @importFrom ggplot2 ggplot geom_histogram labs ggtitle
#'
#' @return a plot of generated data
#'
#' @export
assess.dat.beta <- function(n, q, B, s2, x.thetas, x.gamma, c.gamma,
                            x.correct, c.correct) {

  ## generate data
  dat.list <- gen.data.beta(n = n, q = q, B = B, s2 = s2,
                            x.thetas = x.thetas,
                            x.gamma = x.gamma, c.gamma = c.gamma)

  datf <- dat.list$datf          # full data
  dat0 <- dat.list$dat0          # oracle data
  dat <- dat.list$dat            # observed data
  datcc <- dat.list$datcc        # complete case data

  ## compute empirical values by Z level
  q.hats <- datf %>% group_by(Z) %>% summarise(qhat = mean(X > C)) # cens prop
  x.bars <- datf %>% group_by(Z) %>% summarise(xbar = mean(X))     # mean X

  # estimate distribution of X|Z
  if (x.correct == T) {

    # estimate beta parameters at each level of Z
    x.params.hat <- t(vapply(
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

    # define estimated X|Z density
    eta1 <- function(x, z) {
      dbeta(x = x,
            shape1 = x.params.hat[z + 1, "shape1"],
            shape2 = x.params.hat[z + 1, "shape2"])
    }

  } else {

    # misspecify: estimate marginal beta distribution
    est <- dat %>%
      mutate(left = W,
             right = ifelse(Delta == 1, W, NA)) %>%
      dplyr::select(left, right) %>%
      fitdistrplus::fitdistcens(distr = "beta")
    x.params.hat <- est$estimate

    # define estimated X|Z density
    eta1 <- function(x, z) {
      dbeta(x = x,
            shape1 = x.params.hat["shape1"],
            shape2 = x.params.hat["shape2"])
    }
  }

  # estimate distribution of C|Z
  if (c.correct) {

    # estimate beta parameters at each level of Z
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

    # define estimated C|Z density
    eta2 <- function(c, z) {
      dbeta(x = c,
            shape1 = c.params.hat[z + 1, "shape1"],
            shape2 = c.params.hat[z + 1, "shape2"])
    }

  } else {

    # misspecify: estimate marginal beta distribution
    est <- dat %>%
      mutate(left = W,
             right = ifelse(Delta == 0, W, NA)) %>%
      dplyr::select(left, right) %>%
      fitdistrplus::fitdistcens(distr = "beta")
    c.params.hat <- est$estimate

    # define estimated X|Z density
    eta2 <- function(c, z) {
      dbeta(x = c,
            shape1 = c.params.hat["shape1"],
            shape2 = c.params.hat["shape2"])
    }
  }

  # grid to plot X density
  x.grid <- seq(0, 1, length = 100)

  # data for plotting
  dat <- datf %>%
    mutate(e1 = eta1(x = X, z = Z),
           e2 = eta2(c = C, z = Z))

  # plot data
  plot <- datf %>%
    ggplot() +
    geom_histogram(aes(x = X,
                       y = after_stat(density)),
                   fill = "blue", alpha = 0.5, bins = n/50) +
    geom_histogram(aes(x = C,
                       y = after_stat(density)),
                   fill = "red", alpha = 0.5, bins = n/50) +
    geom_line(aes(x = X,
                  y = e1),
              color = "blue", linewidth = 1) +
    geom_line(aes(x = C,
                  y = e2),
              color = "red", linewidth = 1) +
    facet_wrap(~ Z,
               scales = "free", labeller = label_both) +
    labs(x = "X (blue) or C (red)",
         y = "Density") +
    ggtitle("Estimated vs Observed Distributions of X|Z and C|Z",
            subtitle = paste0("q.hats = ", paste(round(q.hats$qhat, 2), collapse = ", "), "; ",
                              "x.bars = ", paste(round(x.bars$xbar, 2), collapse = ", ")))

  return(list(dat = dat,
              plot = plot))
}

