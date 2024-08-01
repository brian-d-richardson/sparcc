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
                         labels = c("Oracle", "Complete Case", "MLE", "SPARCC")),
         nuis.x = factor(str_sub(str_sub(
           gsub("[.]", "", method.param), 4, -2), 1, 1),
           levels = c("", "1", "0",  "2"),
           labels = c("-", "X Correct", "X Incorrect", "X Nonpar")),
         nuis.c = factor(str_sub(str_sub(
           gsub("[.]", "", method.param), 4, -2), 2, 2),
           levels = c("", "1", "0",  "2"),
           labels = c("-", "C Correct", "C Incorrect", "C Nonpar")),
         param = factor(str_sub(method.param, -1, -1)),
         B.true = B.true[param]) %>%
  filter(name == "B")

# prep for plots ----------------------------------------------------------

param.labs <- c("\u03B20", "\u03B21", "log\u03C3\u00B2")

pal_light <- c('#EE6677', '#228833', '#4477AA', '#CCBB44',
               '#66CCEE','#AA3377', '#BBBBBB')


nuis.x <- unique(sim.out.long$nuis.x)
nuis.c <- unique(sim.out.long$nuis.c)

# create plots ------------------------------------------------------------

# Oracle and Complete Case
seq1 <- ggplot(
  filter(sim.out.long,
         param == 2,
         q == 0.8,
         method %in% c("Oracle", "Complete Case")),
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
        legend.position = "bottom") +
  labs(y = "") +
  guides(fill = "none", color = "none")

# Oracle, Complete Case, and MLE
seq2 <- ggplot(
  filter(sim.out.long,
         param == 2,
         q == 0.8,
         method %in% c("Oracle", "Complete Case", "MLE")),
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
        legend.position = "bottom") +
  labs(y = "") +
  guides(fill = "none", color = "none")

# Complete Case and SPARCC
seq3 <- ggplot(
  filter(sim.out.long,
         param == 2,
         q == 0.8,
         method %in% c("Complete Case", "SPARCC")),
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
        legend.position = "bottom") +
  labs(y = "") +
  guides(fill = "none", color = "none")

# save plots --------------------------------------------------------------

seq1
ggsave("create_figures/robsim_seq1.png",
       dpi = 300,
       width = 6,
       height = 4)

seq2
ggsave("create_figures/robsim_seq2.png",
       dpi = 300,
       width = 6,
       height = 4)

seq3
ggsave("create_figures/robsim_seq3.png",
       dpi = 300,
       width = 6,
       height = 4)
