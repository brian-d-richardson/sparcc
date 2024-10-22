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

# true beta
B.true <- c(1, 10, 2, log(1))

# load simulation results from each of 10 clusters
sim.out.list <- lapply(
  X = 0:9,
  FUN = function(clust) {
    cbind(clust,
          read.csv(paste0(
            "simulation/sim1_data/sd",
            clust, ".csv")))
  })

# combine simulation results into 1 data frame
sim.out <- bind_rows(sim.out.list)

colnames(sim.out) <- gsub(".B", ".", colnames(sim.out))
colnames(sim.out) <- gsub("B2", "trueB2", colnames(sim.out))

# make long data frame
sim.out.long <- sim.out %>%
  pivot_longer(cols = starts_with(c("B", "V")),
               names_to = "method.param",
               values_to = "estimate") %>%
  mutate(estimate = ifelse(abs(estimate) > 100, NA, estimate),
         name = factor(str_sub(method.param, 1, 1)),
         method = factor(str_sub(method.param, 3, 4),
                         levels = c("or", "cc", "ml", "sp"),
                         labels = c("Oracle", "CC", "MLE", "SPARCC")),
         nuis.x = factor(str_sub(str_sub(
           gsub("[.]", "", method.param), 4, -2), 1, 1),
           levels = c("", "1", "0",  "2"),
           labels = c("-", "X Right", "X Wrong", "X Nonpar")),
         nuis.c = factor(str_sub(str_sub(
           gsub("[.]", "", method.param), 4, -2), 2, 2),
           levels = c("", "1", "0",  "2"),
           labels = c("-", "C Right", "C Wrong", "C Nonpar")),
         param = factor(str_sub(method.param, -1, -1)),
         B.true = B.true[param]) %>%
  filter(name == "B")

# prep for plots ----------------------------------------------------------

pal_light <- c('#EE6677', '#228833', '#4477AA', '#CCBB44',
               '#66CCEE','#AA3377', '#BBBBBB')


nuis.x <- unique(sim.out.long$nuis.x)
nuis.c <- unique(sim.out.long$nuis.c)

font.size <- 12

# create plots ------------------------------------------------------------

# Oracle and CC
seq1 <- ggplot(
  filter(sim.out.long,
         param == 2,
         q == 0.8,
         method %in% c("Oracle", "CC")),
  aes(y = estimate,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 0.75) +
  facet_nested(~ method,
               scales = "free") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = font.size),
        strip.text.x = element_text(size = font.size)) +
  labs(y = "Estimated Coefficient") +
  guides(fill = "none", color = "none") +
  scale_fill_manual(values = pal_light[c(1, 2, 3)])

# Oracle, CC, and MLE
seq2 <- ggplot(
  filter(sim.out.long,
         param == 2,
         q == 0.8,
         method %in% c("Oracle", "CC", "MLE")),
  aes(y = estimate,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 0.75) +
  facet_nested(~ method + nuis.x,
               scales = "free") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = font.size),
        strip.text.x = element_text(size = font.size - 2)) +
  labs(y = "Estimated Coefficient") +
  guides(fill = "none", color = "none") +
  scale_fill_manual(values = pal_light[c(1, 2, 3)])

# SPARCC
seq3 <- ggplot(
  filter(sim.out.long,
         param == 2,
         q == 0.8,
         method %in% c("SPARCC")),
  aes(y = estimate,
      fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 0.75) +
  facet_nested(~ method + nuis.x + nuis.c,
               scales = "free") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = font.size),
        strip.text.x = element_text(size = font.size)) +
  labs(y = "Estimated Coefficient") +
  guides(fill = "none", color = "none") +
  scale_fill_manual(values = pal_light[5])

# all estimators
seq4 <- sim.out.long %>%
  #mutate(method_group = ifelse(method %in% levels(method)[1:3], 0, 1)) %>%
  filter(param == 2,
         q == 0.8) %>%
  ggplot(aes(y = estimate,
             fill = method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = B.true),
             linetype = "dashed",
             linewidth = 0.75) +
  facet_nested(~ method + nuis.x + nuis.c,
               scales = "free") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = font.size),
        strip.text.x = element_text(size = font.size)) +
  labs(y = "Estimated Coefficient") +
  guides(fill = "none", color = "none") +
  scale_fill_manual(values = pal_light[c(1, 2, 3, 5)])

# save plots --------------------------------------------------------------

seq1
ggsave("create_figures/robsim_seq1.png",
       dpi = 300,
       width = 6,
       height = 3)

seq2
ggsave("create_figures/robsim_seq2.png",
       dpi = 300,
       width = 6,
       height = 3)

seq3
ggsave("create_figures/robsim_seq3.png",
       dpi = 300,
       width = 6,
       height = 3)

seq4
ggsave("create_figures/robsim_seq4.png",
       dpi = 300,
       width = 9,
       height = 3)


# create plots with bias only ---------------------------------------------

# summarize bias by setting and method
sim.out.bias <- sim.out.long %>%
  filter(name == "B") %>%
  group_by(n, q, method, nuis.x, nuis.c) %>%
  summarise(bias = mean(estimate - B.true, na.rm = T))

# plot all methods
ggplot(filter(sim.out.bias, q == 0.8),
       aes(y = bias,
           color = nuis.x,
           x = method,
           shape = nuis.x)) +
  geom_point()












