###############################################################################
###############################################################################

# Simulation Setup Figures

# Brian Richardson

# 2024-03-04

# Purpose: create figure to illustrate simulation setup

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
setwd("C:/Users/Brian Richardson/OneDrive - University of North Carolina at Chapel Hill/Desktop/Garcia Lab/Research/Projects-Active/HDSP/sparcc")
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(devtools)
load_all()

# define parameters -------------------------------------------------------

n <- 10000            # sample size
q <- 0.8              # censoring proportion
B <- c(1, 2)          # beta
s2 <- 1.1             # variance of Y|X,Z
x.mean <- 1           # mean of X


# create plots ------------------------------------------------------------

# both gamma
plot.11 <- assess.dat(
  n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
  x.shape = 2, c.shape = 2,
  specify.x.gamma = T, specify.c.gamma = T) +
  labs(x = "", y = "") +
  ggtitle(element_blank(), subtitle = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(0, 2.8)

# X gamma, C exp
plot.10 <- assess.dat(
  n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
  x.shape = 2, c.shape = 2,
  specify.x.gamma = T, specify.c.gamma = F) +
  labs(x = "", y = "") +
  ggtitle(element_blank(), subtitle = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(0, 2.8)

# X exp, C gamma
plot.01 <- assess.dat(
  n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
  x.shape = 2, c.shape = 2,
  specify.x.gamma = F, specify.c.gamma = T) +
  labs(x = "", y = "") +
  ggtitle(element_blank(), subtitle = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(0, 2.8)

# both exp
plot.00 <- assess.dat(
  n = n, q = q, B = B, s2 = s2, x.mean = x.mean,
  x.shape = 2, c.shape = 2,
  specify.x.gamma = F, specify.c.gamma = F) +
  labs(x = "", y = "") +
  ggtitle(element_blank(), subtitle = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(0, 2.8)

plotall <- grid.arrange(plot.11, plot.10, plot.01, plot.00)

# save plot ---------------------------------------------------------------

ggsave("create_figures/sim_setup_11.png", plot = plot.11, dpi = 300,
       width = 5, height = 4)
ggsave("create_figures/sim_setup_10.png", plot = plot.10, dpi = 300,
       width = 5, height = 4)
ggsave("create_figures/sim_setup_01.png", plot = plot.01, dpi = 300,
       width = 5, height = 4)
ggsave("create_figures/sim_setup_00.png", plot = plot.00, dpi = 300,
       width = 5, height = 4)

