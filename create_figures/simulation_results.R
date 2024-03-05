###############################################################################
###############################################################################

# Simulation Result Figures

# Brian Richardson

# 2024-02-29

# Purpose: create figures to illustrate simulation results

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
setwd("C:/Users/Brian Richardson/OneDrive - University of North Carolina at Chapel Hill/Desktop/Garcia Lab/Research/Projects-Active/HDSP/sparcc")
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyverse)
library(devtools)
load_all()

# load data ---------------------------------------------------------------

# true (beta, log s2)
B <- c(1, 2, log(0.81))

# load simulation results from each of 10 clusters
sim.out.list <- lapply(
  X = 0:9,
  FUN = function(clust) {
    cbind(clust,
          read.csv(paste0("simulation/sim_data/sim1/unownvar_knownetas/sd",
                          clust, ".csv")))
  })


# combine simulation results into 1 data frame
sim.out <- bind_rows(sim.out.list)

# make long data frame
sim.out.long <- sim.out %>%
  pivot_longer(cols = starts_with("B"),
               names_to = "method.param",
               values_to = "estimate") %>%
  mutate(method = factor(substr(method.param, 2, 3),
                         levels = c("or", "cc", "ml", "sp")),
         param = factor(substr(method.param, 4, 4)),
         B.true = B[param])

# extract simulation parameters
q <- unique(sim.out$q)
n <- unique(sim.out$n)
x.shape <- unique(sim.out$x.shape)
c.shape <- unique(sim.out$c.shape)
mx <- unique(sim.out$mx)
mc <- unique(sim.out$mc)
my <- unique(sim.out$my)
n.rep <- nrow(sim.out) /
  n_distinct(dplyr::select(sim.out, q, n, x.shape, c.shape, mx, mc, my))

# make labels for plots
method.labs <- c("Oracle",
                 "Complete Case",
                 "MLE",
                 "Semiparametric")

names(method.labs) <- c("or", "cc", "ml", "sp")

q.labs <- paste0("q = ", q)
names(q.labs) <- q

n.labs <- paste0("n = ", n)
names(n.labs) <- n

shapex.labs <- c("X Correct", "X Incorrect")
names(shapex.labs) <- x.shape

shapec.labs <- c("C Correct", "C Incorrect")
names(shapec.labs) <- c.shape

mx.labs <- paste0("mx = ", mx)
names(mx.labs) <- mx

mc.labs <- paste0("mc = ", mc)
names(mc.labs) <- mc

my.labs <- paste0("my = ", my)
names(my.labs) <- my

y.range <- sim.out.long %>%
  filter(param == 2,
         q == 0.8,
         x.shape == 1,
         c.shape == 1,
         method == "cc") %>%
  reframe(r = range(estimate)) %>%
  c()

param.labs <- c("\u03B20", "\u03B21", "log\u03C3\u00B2")

# color palettes ----------------------------------------------------------

pal_light <- c('#EE6677', '#228833', '#4477AA', '#CCBB44',
               '#66CCEE','#AA3377', '#BBBBBB')
pal_dark <- c('#991122', '#114419', '#223b55', '#6b611d',
              '#117799', '#55193b', '#5d5d5d')

pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
         "#0072B2", "#D55E00", "#CC79A7")

font.size <- 18

# plot oracle data --------------------------------------------------------

## oracle and complete case only
plot1 <- ggplot(
  filter(sim.out.long,
         param == 2,
         q == 0.8,
         x.shape == 1,
         c.shape == 1,
         method %in% c("or", "cc")),
  aes(y = estimate,
      color = method,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 0.6,
             color = pal_light[4]) +
  labs(y = "",
       fill = "",
       color = "") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = font.size),
        legend.position = "bottom") +
  scale_fill_manual(values = pal_light[c(1, 2, 3, 6)],
                    labels = method.labs) +
  scale_color_manual(values = pal_dark[c(1, 2, 3, 6)],
                     labels = method.labs) +
  ylim(y.range$r)
plot1

## oracle, complete case, and MLE
plot2 <- ggplot(
  filter(sim.out.long,
         param == 2,
         q == 0.8,
         x.shape == 1,
         c.shape == 1,
         method %in% c("or", "cc", "ml")),
  aes(y = estimate,
      color = method,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 0.6,
             color = pal_light[4]) +
  labs(y = "",
       fill = "",
       color = "") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = font.size),
        legend.position = "bottom") +
  scale_fill_manual(values = pal_light[c(1, 2, 3, 6)],
                    labels = method.labs) +
  scale_color_manual(values = pal_dark[c(1, 2, 3, 6)],
                     labels = method.labs) +
  ylim(y.range$r)
plot2

## oracle, complete case, and MLE with misspecified X
plot3 <- ggplot(
  filter(sim.out.long,
         param == 2,
         q == 0.8,
         c.shape == 1,
         method %in% c("or", "cc", "ml")),
  aes(y = estimate,
      color = method,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 0.6,
             color = pal_light[4]) +
  facet_grid(~ x.shape,
             scales = "free",
             labeller = labeller(q = q.labs,
                                 x.shape = shapex.labs,
                                 c.shape = shapec.labs)) +
  labs(y = "",
       fill = "",
       color = "") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = font.size),
        legend.position = "bottom",
        strip.text.x = element_text(size = font.size)) +
  scale_fill_manual(values = pal_light[c(1, 2, 3, 6)],
                    labels = method.labs) +
  scale_color_manual(values = pal_dark[c(1, 2, 3, 6)],
                     labels = method.labs) +
  ylim(y.range$r)
plot3

## oracle, complete case, MLE, and SP with misspecified X
plot4 <- ggplot(
  filter(sim.out.long,
         param == 2,
         q == 0.8,
         c.shape == 1),
  aes(y = estimate,
      color = method,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 0.6,
             color = pal_light[4]) +
  facet_grid(~ x.shape,
             scales = "free",
             labeller = labeller(q = q.labs,
                                 x.shape = shapex.labs,
                                 c.shape = shapec.labs)) +
  labs(y = "",
       fill = "",
       color = "") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = font.size),
        legend.position = "bottom",
        strip.text.x = element_text(size = font.size)) +
  scale_fill_manual(values = pal_light[c(1, 2, 3, 6)],
                    labels = method.labs) +
  scale_color_manual(values = pal_dark[c(1, 2, 3, 6)],
                     labels = method.labs) +
  ylim(y.range$r)
plot4

## oracle, complete case, MLE, and SP with misspecified X, C
plot5 <- ggplot(
  filter(sim.out.long,
         param == 2,
         q == 0.8),
  aes(y = estimate,
      color = method,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 0.6,
             color = pal_light[4]) +
  facet_nested(~ x.shape + c.shape,
               scales = "free",
               labeller = labeller(q = q.labs,
                                   x.shape = shapex.labs,
                                   c.shape = shapec.labs)) +
  labs(y = "",
       fill = "",
       color = "") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = font.size),
        legend.position = "bottom",
        strip.text.x = element_text(size = font.size)) +
  scale_fill_manual(values = pal_light[c(1, 2, 3, 6)],
                    labels = method.labs) +
  scale_color_manual(values = pal_dark[c(1, 2, 3, 6)],
                     labels = method.labs) +
  ylim(y.range$r)
plot5

# save plots
ggsave("create_figures/sim_1.png", plot = plot1, dpi = 300,
       width = 7, height = 4)
ggsave("create_figures/sim_2.png", plot = plot2, dpi = 300,
       width = 7, height = 4)
ggsave("create_figures/sim_3.png", plot = plot3, dpi = 300,
       width = 7, height = 4)
ggsave("create_figures/sim_4.png", plot = plot4, dpi = 300,
       width = 7, height = 4)
ggsave("create_figures/sim_5.png", plot = plot5, dpi = 300,
       width = 7, height = 4)

