###############################################################################
###############################################################################

# Fake ENROLL-HD Data Analysis Results for All Outcomes

# Anonymous

# 2025-08-26

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(kableExtra)
library(here)
setwd(here())

# load data ---------------------------------------------------------------

res.tms <- read.csv("derived-data/tms-results/res.csv")
res.sdmt <- read.csv("derived-data/sdmt-results/res.csv")
res.cuhdrs <- read.csv("derived-data/cUHDRS-results/res.csv")
dat <- read.csv("derived-data/dat_cuhdrs.csv")

res <- rbind(cbind(outcome = "TMS", res.tms),
             cbind(outcome = "SDMT", res.sdmt),
             cbind(outcome = "cUHDRS", res.cuhdrs)) %>%
  mutate(Method = ifelse(method == "Semipar. (Param.)", "Semipar",
                  ifelse(method == "Semipar. (Nonpar.)", "Semipar",
                  ifelse(method == "MLE", "MLE", "CC"))),
         Method = factor(Method, levels = c("Semipar", "MLE", "CC")),
         nuisX = ifelse(method %in% c("Semipar. (Param.)", "MLE"), "Param",
                 ifelse(method == "Semipar. (Nonpar.)", "Nonpar", "-")),
         nuisX = factor(nuisX, levels = c("Param", "Nonpar", "-")),
         nuisC = ifelse(method == "Semipar. (Param.)", "Param",
                 ifelse(method == "Semipar. (Nonpar.)", "Nonpar", "-")),
         nuisC = factor(nuisC, levels = c("Param", "Nonpar", "-")),
         outcome = factor(outcome, levels = c("TMS", "SDMT", "cUHDRS")))

# prep for plots ----------------------------------------------------------

# underscore (\u2081) not working in plot
param.labs <- c("\u03B21", "\u03B22", "\u03B23", "\u03B24", "log\u03C3\u00B2")
names(param.labs) <- 1:5

method.labs <- c("SPARCC (Param)", "SPARCC (Nonpar)", "MLE", "CC")
names(method.labs) <- c("sp_param", "sp_nonpar", "ml", "cc")

# view results ------------------------------------------------------------

## Plot estimates and CIs
ggplot(
  data = res,
  aes(x = method,
      y = est,
      ymin = ci.lower,
      ymax = ci.upper,
      color = method)) +
  geom_point() +
  geom_errorbar() +
  facet_grid(outcome ~ param,
             labeller = labeller(param = param.labs),
             scales = "free_y") +
  labs(y = "Estimate with 95% CI",
       color = "Method") +
  ggtitle("ENROLL-HD Analysis Results",
          subtitle = "Outcome ~ (Time to Diagnosis) + CAP Category") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

## Table of results
tbl1 <- res %>%
  mutate(ci = paste0("(", round(ci.lower, 2),
                     ", ", round(ci.upper, 2), ")")) %>%
  dplyr::select(outcome, param, Method, nuisX, nuisC, est, ste, ci) %>%
  arrange(outcome, param, Method, nuisX) %>%
  mutate_if(is.numeric, function(x) round(x, 2))

tbl1 %>%
  kable(format = "latex",
        digits = 2,
        align = c(rep("c", 4),
                  rep("r", 3)),
        booktabs = TRUE,
        linesep = c("", "", "", "\\addlinespace\\hline"),
        escape = FALSE,
        col.names = c("Outcome", "Parameter", "Method", "X|Z", "C|Z",
                      "Estimate", "Std Error", "95% CI")) %>%
  kable_styling("striped") %>%
  row_spec(row = 0, bold = TRUE)

## Table of results with all params and with est(se) only
tbl2 <- res %>%
  mutate(sig = ifelse(abs(est / ste) > qnorm(0.975), "*", ""),
         est.ste = paste0(formatC(est, format = "f", digits = 2), "(",
                          formatC(ste, format = "f", digits = 2), ")",
                          sig)) %>%
  dplyr::select(outcome, Method, nuisX, nuisC, param, est.ste) %>%
  pivot_wider(values_from = est.ste,
              names_from = param,
              id_cols = c(outcome, Method, nuisX, nuisC))

tbl2 %>%
  kable(format = "latex",
        digits = 2,
        align = c(rep("c", 4),
                  rep("r", 5)),
        booktabs = TRUE,
        linesep = c("", "", "", "\\addlinespace\\hline"),
        escape = FALSE,
        col.names = c("Outcome", "Estimator", "X|Z", "C|Z",
                      "beta_1", "beta_2", "beta_3",
                      "beta_4", "log(sigma2)")) %>%
  add_header_above(c(" " = 4,
                     "Est(Std Error)" = 5)) %>%
  kable_styling("striped") %>%
  row_spec(row = 0, bold = TRUE)

# plot regression lines ---------------------------------------------------

## get predicted values over range of X values
len.new <- 100
Xnew <- data.frame(
  int = 1,
  X = seq(min(dat$W), max(dat$W), length = len.new),
  Z = rep(0:1, each = len.new)) %>%
  mutate(XZ = X * Z) %>%
  as.matrix()

Yhats <- do.call(rbind, lapply(
  unique(res$outcome),
  function(outcome) {
    Yhat.list <- lapply(
      names(method.labs),
      function(m) {

        BV <- readRDS(paste0("derived-data/", outcome,
                             "-results/", m, "_res"))
        B <- BV$B[1:4]; V <- BV$V[1:4, 1:4]
        g <- ifelse(outcome == "TMS",
                    function(x) exp(x) - 1,
                    function(x) x)

        Yhat <- Xnew %*% B
        Yhat.se <- sqrt(diag(Xnew %*% V %*% t(Xnew)))
        pred <- data.frame(
          method = m,
          int = Xnew[,1],
          X = Xnew[,2],
          Z = factor(Xnew[,3]),
          Yhat = g(Yhat),
          Ylower = g(Yhat - qnorm(0.975) * Yhat.se),
          Yupper = g(Yhat + qnorm(0.975) * Yhat.se)
        )

      })
    return(cbind(outcome = outcome, do.call(rbind, Yhat.list)))
  })) %>%
  mutate(Xscaled = (X * (2796.01 - 17.99) + 17.99) / 365,
         Z = factor(Z,
                    levels = 0:1,
                    labels = c("Low-Medium", "High")),
         method = factor(method,
                         levels = c("sp_param", "sp_nonpar", "ml", "cc")))

## plot fitted lines and 95% confidence intervals
ggplot(data = Yhats,
       aes(x = Xscaled,
           y = Yhat,
           ymin = Ylower,
           ymax = Yupper,
           color = Z,
           fill = Z)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  facet_grid(outcome ~ method,
             scales = "free",
             labeller = labeller(method = method.labs)) +
  theme_bw() +
  scale_x_reverse() +
  labs(y = "Estimated Mean Outcome (and 95% Confidence Interval)",
       x = "Years Until Diagnosis",
       color = "CAP Score:",
       fill = "CAP Score:") +
  theme(legend.position = "bottom")

ggsave("figures/analysis_fig1.png", dpi = 300, width = 6, height = 6)


## plot fitted lines and 95% confidence intervals for cUHDRS only
Yhats %>%
  filter(outcome == "cUHDRS") %>%
  mutate(method = factor(
    method,
    levels = c("cc", "ml", "sp_param", "sp_nonpar"),
    labels = c("Complete Case", "MLE", "SPARCC (Param)", "SPARCC (Nonpar)"))) %>%
  ggplot(
    aes(x = Xscaled,
        y = Yhat,
        ymin = Ylower,
        ymax = Yupper,
        color = Z,
        fill = Z)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  facet_grid(~ method,
             scales = "free") +
  theme_bw() +
  scale_x_reverse() +
  labs(y = "cUHDRS",
       x = "Years Until Diagnosis",
       color = "CAP Score:",
       fill = "CAP Score:") +
  theme(legend.position = "bottom")

ggsave("figures/analysis_fig2.png", dpi = 300, width = 6, height = 3)
