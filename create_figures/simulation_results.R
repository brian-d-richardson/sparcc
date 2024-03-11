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
B <- c(1, 0.2, log(0.09))

# load simulation results from each of 10 clusters
sim.out.list <- lapply(
  X = 0:9,
  FUN = function(clust) {
    cbind(clust,
          read.csv(paste0("simulation/sim_data/sim1/fine_search_4/sd",
                          clust, ".csv")))
  })


# combine simulation results into 1 data frame
sim.out <- bind_rows(sim.out.list)
colnames(sim.out) <- gsub(".B", "", colnames(sim.out))

# make long data frame
sim.out.long <- sim.out %>%
  select(!B2) %>%
  pivot_longer(cols = starts_with("B"),
               names_to = "method.param",
               values_to = "estimate") %>%
  mutate(method = factor(substr(method.param, 2, 3),
                         levels = c("or", "cc", "ml", "sp")),
         param = factor(substr(method.param, 4, 4)),
         specify.x.gamma = factor(specify.x.gamma,
                                  levels = c(1, 0)),
         specify.c.gamma = factor(specify.c.gamma,
                                  levels = c(1, 0)),
         B.true = B[param])

# extract simulation parameters
q <- unique(sim.out$q)
n <- unique(sim.out$n)
specify.x.gamma <- unique(sim.out$specify.x.gamma)
specify.c.gamma <- unique(sim.out$specify.c.gamma)
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

specx.labs <- c("X Correct", "X Incorrect")
names(specx.labs) <- specify.x.gamma

specc.labs <- c("C Correct", "C Incorrect")
names(specc.labs) <- specify.c.gamma

mx.labs <- paste0("mx = ", mx)
names(mx.labs) <- mx

mc.labs <- paste0("mc = ", mc)
names(mc.labs) <- mc

my.labs <- paste0("my = ", my)
names(my.labs) <- my

y.range <- sim.out.long %>%
  filter(param == 2,
         q == 0.7,
         specify.x.gamma == 1,
         specify.c.gamma == 1,
         method == "cc") %>%
  reframe(r = range(estimate)) %>%
  c()

param.labs <- c("\u03B20", "\u03B21", "log\u03C3\u00B2")


# summary data for annotation ---------------------------------------------

sum.dat <- sim.out.long %>%
  filter(method != "ml" |                 # filter extreme MLE observations
         abs(estimate - B.true) < 1) %>%
  group_by(n, q, specify.x.gamma, specify.c.gamma, method, param) %>%
  summarize(var = round(1000*var(estimate), 2)) %>%
  mutate(lab = paste("Var = ", var),
         ht = ifelse(method %in% c("or", "ml"), 0.3, 0.1)) %>%
  select(n, q, specify.x.gamma, specify.c.gamma, method, param, lab, ht)

plot.dat <- inner_join(sim.out.long, sum.dat,
                       by = NULL)

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
  filter(plot.dat,
         param == 2,
         q == 0.7,
         specify.x.gamma == 1,
         specify.c.gamma == 1,
         method %in% c("or", "cc")),
  aes(y = estimate,
      x = method,
      color = method,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 1,
             color = pal_light[4]) +
  facet_nested(~ specify.x.gamma,
               scales = "free",
               labeller = labeller(n = n.labs,
                                   specify.x.gamma = specx.labs,
                                   specify.c.gamma = specc.labs)) +
  labs(x = "",
       y = "",
       fill = "",
       color = "") +
  geom_text(aes(x = method, label = lab, y = ht)) +
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
plot1

## oracle, complete case, and MLE
plot2 <- ggplot(
  filter(plot.dat,
         param == 2,
         q == 0.7,
         specify.x.gamma == 1,
         specify.c.gamma == 1,
         method %in% c("or", "cc", "ml")),
  aes(y = estimate,
      x = method,
      color = method,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 1,
             color = pal_light[4]) +
  facet_nested(~ specify.x.gamma,
               scales = "free",
               labeller = labeller(n = n.labs,
                                   specify.x.gamma = specx.labs,
                                   specify.c.gamma = specc.labs)) +
  labs(x = "",
       y = "",
       fill = "",
       color = "") +
  geom_text(aes(x = method, label = lab,
                y = ht)) +
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
plot2

## oracle, complete case, and MLE with misspecified X
plot3 <- ggplot(
  filter(plot.dat,
         param == 2,
         q == 0.7,
         specify.c.gamma == 1,
         method %in% c("or", "cc", "ml")),
  aes(y = estimate,
      x = method,
      color = method,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 1,
             color = pal_light[4]) +
  facet_nested(~ specify.x.gamma,
               scales = "free",
               labeller = labeller(n = n.labs,
                                   specify.x.gamma = specx.labs,
                                   specify.c.gamma = specc.labs)) +
  labs(x = "",
       y = "",
       fill = "",
       color = "") +
  geom_text(aes(x = method, label = lab, y = ht)) +
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
  filter(plot.dat,
         param == 2,
         q == 0.7,
         specify.c.gamma == 1),
  aes(y = estimate,
      x = method,
      color = method,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 1,
             color = pal_light[4]) +
  facet_nested(~ specify.x.gamma,
               scales = "free",
               labeller = labeller(n = n.labs,
                                   specify.x.gamma = specx.labs,
                                   specify.c.gamma = specc.labs)) +
  labs(x = "",
       y = "",
       fill = "",
       color = "") +
  geom_text(aes(x = method, label = lab, y = ht)) +
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
  filter(plot.dat,
         param == 2,
         q == 0.7),
  aes(y = estimate,
      color = method,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 0.75,
             color = pal_light[4]) +
  facet_nested(~ specify.x.gamma + specify.c.gamma,
               scales = "free",
               labeller = labeller(n = n.labs,
                                   specify.x.gamma = specx.labs,
                                   specify.c.gamma = specc.labs)) +
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

