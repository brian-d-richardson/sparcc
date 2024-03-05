###############################################################################
###############################################################################

# Basic Censoring Figures

# Brian Richardson

# 2024-02-27

# Purpose: create figure to illustrate basic censoring

###############################################################################
###############################################################################

# prep workspace ----------------------------------------------------------

rm(list = ls())
#setwd(dirname(getwd()))
setwd("C:/Users/Brian Richardson/OneDrive - University of North Carolina at Chapel Hill/Desktop/Garcia Lab/Research/Projects-Active/HDSP/sparcc")
library(ggplot2)
library(dplyr)
library(tidyverse)
library(devtools)
load_all()

# define parameters -------------------------------------------------------

n <- 200              # sample size
q <- 0.8              # censoring proportion
B <- c(1, 2)          # beta
s2 <- 1.1             # variance of Y|X,Z
x.mean <- 1           # mean of X
set.seed(2)

# generate data -----------------------------------------------------------

x.shape <- 4
c.shape <- 4
x.rate <- x.shape / x.mean # rate parameter for gamma distribution of X
c.rate <- get.c.rate(      # rate parameter for gamma distribution of C
  q = q,
  x.rate = x.rate,
  x.shape = x.shape,
  c.shape = c.shape)

dat.list <- gen.data(n = n, q = q, B = B, s2 = s2,
                     x.rate = x.rate, x.shape = x.shape,
                     c.rate = c.rate, c.shape = c.shape)

datf <- dat.list$datf          # full data
dat0 <- dat.list$dat0          # oracle data
dat <- dat.list$dat            # observed data
datcc <- dat.list$datcc        # complete case data


# add fitted values -------------------------------------------------------

lm.oracle <- lm(Y ~ X, data = datf) # oracle linear model
lm.naive <- lm(Y ~ W, data = dat)   # naive linear model
lm.cc <- lm(Y ~ W, data = datcc)    # complete case
true.fit <- B[1] + B[2]*datf$X

plot.dat <- cbind(
  datf,
  select(dat, W, Delta),
  yhat.oracle = predict(lm.oracle, datf,
                               interval = "confidence"),
  yhat.naive = predict(lm.naive, rename(datf, W = "X"),
                              interval = "confidence"),
  yhat.cc = predict(lm.cc, rename(datf, W = "X"),
                           interval = "confidence")) %>%
  pivot_longer(cols = starts_with("yhat"),
               names_to = "method.val",
               values_to = "pred") %>%
  mutate(method.val = gsub("yhat.", "", method.val)) %>%
  mutate(method = str_split_i(method.val, "[.]", 1),
         val = str_split_i(method.val, "[.]", 2)) %>%
  select(!method.val) %>%
  pivot_wider(names_from = val,
              values_from = pred)

# color palettes ----------------------------------------------------------

pal_light <- c('#EE6677', '#228833', '#4477AA', '#CCBB44',
               '#66CCEE','#AA3377', '#BBBBBB')
pal_dark <- c('#991122', '#114419', '#223b55', '#6b611d',
              '#117799', '#55193b', '#5d5d5d')

pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
         "#0072B2", "#D55E00", "#CC79A7")

font.size <- 16
pt.size <- 3
line.width <- 1.5

# plot oracle data --------------------------------------------------------

# plot true regression line
plot1 <- plot.dat %>%
  filter(method == "oracle") %>%
  ggplot(aes(x = X, y = Y)) +
  geom_point(color = pal[3],
             size = pt.size) +
  geom_line(aes(y = true.fit),
            color = "black",
            linewidth = line.width,
            linetype = "dashed") +
  xlim(range(datf$X)) +
  ylim(range(datf$Y)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size))
plot1

# add oracle regression line
plot2 <- plot.dat %>%
  filter(method == "oracle") %>%
  ggplot(aes(x = X, y = Y)) +
  geom_point(color = pal[3],
             size = pt.size) +
  geom_line(aes(y = true.fit),
            color = "black",
            linewidth = line.width,
            linetype = "dashed") +
  geom_line(aes(y = fit),
            color = pal[8],
            linewidth = line.width) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              fill = pal[8],
              alpha = 0.3) +
  xlim(range(datf$X)) +
  ylim(range(datf$Y)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size))
plot2

# indicate censored observations
plot3 <- plot.dat %>%
  filter(method == "oracle") %>%
  ggplot(aes(x = X, y = Y)) +
  geom_point(aes(x = W),
             color = pal[3],
             size = pt.size) +
  geom_point(data = filter(plot.dat, method == "oracle", Delta == 0),
             shape = 13,
             color = "red",
             size = pt.size) +
  geom_line(aes(y = true.fit),
            color = "black",
            linewidth = line.width,
            linetype = "dashed") +
  xlim(range(datf$X)) +
  ylim(range(datf$Y)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size))
plot3

# add naive regression line
plot4 <- plot.dat %>%
  filter(method == "naive") %>%
  ggplot(aes(x = X, y = Y)) +
  geom_point(aes(x = W),
             color = pal[3],
             size = pt.size) +
  geom_line(aes(y = true.fit),
            color = "black",
            linewidth = line.width,
            linetype = "dashed") +
  geom_line(aes(y = fit),
            color = pal[8],
            linewidth = line.width) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              fill = pal[8],
              alpha = 0.3) +
  xlim(c(max(min(datf$X), min(dat$W)),
         min(max(datf$X), max(dat$W)))) +
  ylim(range(datf$Y)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size))
plot4

# limit to complete cases regression line
plot5 <- plot.dat %>%
  filter(method == "cc") %>%
  ggplot(aes(x = X, y = Y)) +
  geom_point(data = filter(plot.dat, method == "oracle", Delta == 1),
             aes(x = W),
             color = pal[3],
             size = pt.size) +
  geom_line(aes(y = true.fit),
            color = "black",
            linewidth = line.width,
            linetype = "dashed") +
  geom_line(aes(y = fit),
            color = pal[8],
            linewidth = line.width) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              fill = pal[8],
              alpha = 0.2) +
  xlim(range(datcc$W)) +
  ylim(range(datf$Y)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size))
plot5

# fixed censoring
q2 <- quantile(datf$X, 0.2)
plot6 <- plot.dat %>%
  filter(method == "oracle") %>%
  ggplot(aes(x = X, y = Y)) +
  geom_rect(aes(xmin = q2, xmax = Inf,
                ymin = -Inf, ymax = Inf),
            fill = "pink", alpha = 0.05) +
  geom_point(data = filter(plot.dat, method == "oracle", X > q2),
             shape = 13,
             color = "red",
             size = pt.size) +
  geom_point(data = filter(plot.dat, method == "oracle", X <= q2),
             aes(x = X),
             color = pal[3],
             size = pt.size) +
  geom_line(aes(y = true.fit),
            color = "black",
            linewidth = line.width,
            linetype = "dashed") +
  geom_vline(xintercept = q2,
             color = "red",
             linewidth = line.width) +
  xlim(range(datf$X)) +
  ylim(range(datf$Y)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size))
plot6

# save plots
ggsave("create_figures/basic_censoring_1.png", plot = plot1, dpi = 300,
       width = 5, height = 4)
ggsave("create_figures/basic_censoring_2.png", plot = plot2, dpi = 300,
       width = 5, height = 4)
ggsave("create_figures/basic_censoring_3.png", plot = plot3, dpi = 300,
       width = 5, height = 4)
ggsave("create_figures/basic_censoring_4.png", plot = plot4, dpi = 300,
       width = 5, height = 4)
ggsave("create_figures/basic_censoring_5.png", plot = plot5, dpi = 300,
       width = 5, height = 4)
ggsave("create_figures/basic_censoring_6.png", plot = plot6, dpi = 300,
       width = 5, height = 4)

