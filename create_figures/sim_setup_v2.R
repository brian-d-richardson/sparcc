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

# sample size and censoring proportion
n <- 8000
q <- 0.8

# outcome model parameters
B <- c(1, 10, 2)
s2 <- 1

# nuisance model parameters
x.thetas <- 0.5 * c(-1, 1)
x.gamma <- 1
c.gamma <- 2

# create plots ------------------------------------------------------------

# X correctly specified
a <- assess.dat.beta(
  n = n, q = q, B = B, s2 = s2,
  x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma,
  x.correct = T, c.correct = T)
adat <- a$dat %>%
  mutate(Z = factor(Z))

aplot <- ggplot(adat,
       aes()) +
  geom_histogram(aes(y = after_stat(density),
                     x = X,
                     fill = Z),
                 bins = 40,
                 alpha = 0.5,
                 position = "identity") +
  geom_line(aes(x = X,
                y = e1,
                color = Z),
            linewidth = 1.5) +
  labs(x = "",
       y = "",
       fill = "Z") +
  scale_fill_manual(values = c("blue", "red")) +
    scale_color_manual(values = c("blue", "red")) +
  theme_bw() +
  guides(color = "none", fill = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

aplot

# X incorrectly specified
b <- assess.dat.beta(
  n = n, q = q, B = B, s2 = s2,
  x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma,
  x.correct = F, c.correct = T)
bdat <- b$dat %>%
  mutate(Z = factor(Z))

bplot <- ggplot(bdat,
       aes()) +
  geom_histogram(aes(y = after_stat(density),
                     x = X,
                     fill = Z),
                 bins = 40,
                 alpha = 0.5,
                 position = "identity") +
  geom_line(aes(x = X,
                y = e1),
            linewidth = 1.5,
            color = "black") +
  labs(x = "",
       y = "",
       fill = "Z") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw() +
  guides(color = "none", fill = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

bplot

# X modeled nonparametrically with B-splines
dat.list <- gen.data.beta(
  n = n, q = q, B = B, s2 = s2,
  x.thetas = x.thetas, x.gamma = x.gamma, c.gamma = c.gamma)
dat <- dat.list$dat
spline.res.x <- fit.spline(dat = dat, m.knots = 5, deg = 3)
eta1.hat <- spline.res.x$dens
cdat <- dat.list$datf %>%
  mutate(eta1.hat = eta1.hat(x = X, z = Z),
         Z = factor(Z))

cplot <- ggplot(cdat, aes()) +
  geom_histogram(aes(y = after_stat(density),
                     x = X,
                     fill = Z),
                 bins = 40,
                 alpha = 0.5,
                 position = "identity") +
  geom_line(aes(x = X,
                y = eta1.hat,
                color = Z),
            linewidth = 1.5) +
  labs(x = "",
       y = "",
       fill = "Z") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw() +
  guides(color = "none", fill = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

cplot

# save plot ---------------------------------------------------------------

ggsave("create_figures/sim_setup_correct.png", plot = aplot, dpi = 300,
       width = 8, height = 4)
ggsave("create_figures/sim_setup_incorrect.png", plot = bplot, dpi = 300,
       width = 8, height = 4)
ggsave("create_figures/sim_setup_nonpar.png", plot = cplot, dpi = 300,
       width = 8, height = 4)

