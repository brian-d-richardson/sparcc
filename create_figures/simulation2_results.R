###############################################################################
###############################################################################

# Simulation Result Figures

# Brian Richardson

# 2024-02-29

# Purpose: create figures to illustrate simulation 2 results

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
setwd("C:/Users/Brian Richardson/OneDrive - University of North Carolina at Chapel Hill/Desktop/Garcia Lab/Research/Projects-Active/HDSP/sparcc")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggh4x)
library(scales)
library(ggthemes)
library(kableExtra)


# load data ---------------------------------------------------------------

# load simulation results from each of 10 clusters
sim.out.list <- lapply(
  X = 0:9,
  FUN = function(clust) {
    cbind(clust,
          read.csv(paste0(
            "simulation/sim2_data/sd",
            clust, ".csv")))
  })

# combine simulation results into 1 data frame
sim.out <- bind_rows(sim.out.list)

colnames(sim.out) <- gsub(".B", "", colnames(sim.out))
colnames(sim.out) <- gsub("B2", "trueB2", colnames(sim.out))

# true beta
B <- c(1, 10, 2, log(1))

# make long data frame
sim.out.long <- sim.out %>%
  pivot_longer(cols = starts_with("B"),
               names_to = "method.param",
               values_to = "estimate") %>%
  mutate(method = factor(substr(method.param, 2, 3),
                         levels = c("or", "cc", "ml", "sp"),
                         labels = c("Oracle",
                                    "Complete Case",
                                    "MLE",
                                    "SPARCC")),
         param = factor(substr(method.param, 4, 4)),
         B.true = B[param])

# prep for plots ----------------------------------------------------------

# extract simulation parameters
q <- unique(sim.out$q)
n.rep <- nrow(sim.out) /
  n_distinct(dplyr::select(sim.out, q))

# make labels for plots
param.labs <- c("\u03B20", "\u03B21", "\u03B22", "log\u03C3\u00B2")
names(param.labs) <- 1:4
font.size <- 16

# colorblind friendly pallette
pal_light <- c('#EE6677', '#228833', '#4477AA', '#CCBB44',
               '#66CCEE', '#AA3377', '#BBBBBB')

# obtain empirical variances by method, parameter, and q
sim.out.sum <- sim.out.long %>%
  group_by(q, method, param) %>%
  filter(!is.na(estimate)) %>%
  summarise(var = var(estimate),
            Var = mean((estimate - B.true) ^ 2),
            se.Var = sd((estimate - B.true) ^ 2) / sqrt(n.rep),
            Var.lower = Var - qnorm(0.975) * se.Var ,
            Var.upper = Var + qnorm(0.975) * se.Var)

# create plots ------------------------------------------------------------

# plot variance over q for beta_1 only, one estimator at a time
seq1 <- ggplot(
  filter(sim.out.sum, param == 2, q <= 0.85,
         method %in% c("Oracle", "Complete Case")),
  aes(x = q,
      y = var,
      ymin = Var.lower,
      ymax = Var.upper,
      fill = method,
      color = method,
      shape = method)) +
  geom_point() +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  labs(x = "Censoring Proportion (q)",
       y = "Variance",
       color = "Method",
       shape = "Method",
       fill = "Method") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = font.size),
        axis.title.x = element_text(size = font.size),
        axis.text.x = element_text(size = font.size - 2),
        legend.text = element_text(size = font.size - 2),
        legend.title = element_text(size = font.size)) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
  scale_color_manual(values = pal_light[c(1, 2, 3, 5)]) +
  scale_fill_manual(values = pal_light[c(1, 2, 3, 5)])

seq2 <- ggplot(
  filter(sim.out.sum, param == 2, q <= 0.85,
         method %in% c("Oracle", "Complete Case", "MLE")),
  aes(x = q,
      y = var,
      ymin = Var.lower,
      ymax = Var.upper,
      fill = method,
      color = method,
      shape = method)) +
  geom_point() +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  labs(x = "Censoring Proportion (q)",
       y = "Variance",
       color = "Method",
       shape = "Method",
       fill = "Method") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = font.size),
        axis.title.x = element_text(size = font.size),
        axis.text.x = element_text(size = font.size - 2),
        legend.text = element_text(size = font.size - 2),
        legend.title = element_text(size = font.size)) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
  scale_color_manual(values = pal_light[c(1, 2, 3, 5)]) +
  scale_fill_manual(values = pal_light[c(1, 2, 3, 5)])

seq3 <- ggplot(
  filter(sim.out.sum, param == 2, q <= 0.85),
  aes(x = q,
      y = var,
      ymin = Var.lower,
      ymax = Var.upper,
      fill = method,
      color = method,
      shape = method)) +
  geom_point() +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  labs(x = "Censoring Proportion (q)",
       y = "Variance",
       color = "Method",
       shape = "Method",
       fill = "Method") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = font.size),
        axis.title.x = element_text(size = font.size),
        axis.text.x = element_text(size = font.size - 2),
        legend.text = element_text(size = font.size - 2),
        legend.title = element_text(size = font.size)) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
  scale_color_manual(values = pal_light[c(1, 2, 3, 5)]) +
  scale_fill_manual(values = pal_light[c(1, 2, 3, 5)])


# save plots --------------------------------------------------------------

seq1
ggsave("create_figures/varsim_seq1.png",
       dpi = 300,
       width = 6,
       height = 4)

seq2
ggsave("create_figures/varsim_seq2.png",
       dpi = 300,
       width = 6,
       height = 4)

seq3
ggsave("create_figures/varsim_seq3.png",
       dpi = 300,
       width = 6,
       height = 4)


