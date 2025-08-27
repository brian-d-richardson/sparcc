###############################################################################
###############################################################################

# ENROLL-HD Fake Data Generation

# Brian Richardson

# 2025-08-26

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

rm(list = ls())
library(dplyr)
library(here)
library(devtools)
setwd(here())
load_all(dirname(getwd()))

# parameters --------------------------------------------------------------

# parameters to emulate in fake date (based on real ENROLL-HD analysis)
n <- 4530 # sample size
q <- 0.819 # censoring proportion

# outcome model coefficients
beta.TMS <- c(1.66, -0.47, 0.94, -0.30)
beta.SDMT <- c(42.22, 11.35, -11.13, 4.94)
beta.cUHDRS <- c(14.27, 2.96, -2.34, 2.02)

# nuisance model coefficients
x.thetas <- 0.5 * c(-1, 1)
x.gamma <- 1
c.gamma <- 2

set.seed(1)

# generate TMS data -------------------------------------------------------

# generate data
data.TMS <- gen.data.beta(
  n = n,
  q = q,
  B = beta.TMS[1:3],
  s2 = exp(beta.TMS[4]),
  x.thetas = x.thetas,
  x.gamma = x.gamma,
  c.gamma = c.gamma)

# select observed data
dat.TMS <- data.TMS$dat

# select full data (with censored X values set to NA)
datf.TMS <- data.TMS$datf %>%
  mutate(Delta = ifelse(X <= C, 1, 0),
         W = ifelse(Delta == 0, C, X),
         X = ifelse(Delta == 0, NA, X))

# save data
write.csv(dat.TMS, "derived-data/dat_tms.csv", row.names = F)
write.csv(datf.TMS, "derived-data/datf_tms.csv", row.names = F)

# generate SDMT data ------------------------------------------------------

data.SDMT <- gen.data.beta(
  n = n,
  q = q,
  B = beta.SDMT[1:3],
  s2 = exp(beta.SDMT[4]),
  x.thetas = x.thetas,
  x.gamma = x.gamma,
  c.gamma = c.gamma)

# select observed data
dat.SDMT <- data.SDMT$dat

# select full data (with censored X values set to NA)
datf.SDMT <- data.SDMT$datf %>%
  mutate(Delta = ifelse(X <= C, 1, 0),
         W = ifelse(Delta == 0, C, X),
         X = ifelse(Delta == 0, NA, X))

# save data
write.csv(dat.SDMT, "derived-data/dat_sdmt.csv", row.names = F)
write.csv(datf.SDMT, "derived-data/datf_sdmt.csv", row.names = F)


# generate cUHDRS data ----------------------------------------------------

data.cUHDRS <- gen.data.beta(
  n = n,
  q = q,
  B = beta.cUHDRS[1:3],
  s2 = exp(beta.cUHDRS[4]),
  x.thetas = x.thetas,
  x.gamma = x.gamma,
  c.gamma = c.gamma)

# select observed data
dat.cUHDRS <- data.cUHDRS$dat

# select full data (with censored X values set to NA)
datf.cUHDRS <- data.cUHDRS$datf %>%
  mutate(Delta = ifelse(X <= C, 1, 0),
         W = ifelse(Delta == 0, C, X),
         X = ifelse(Delta == 0, NA, X))

# save data
write.csv(dat.cUHDRS, "derived-data/dat_cuhdrs.csv", row.names = F)
write.csv(datf.cUHDRS, "derived-data/datf_cuhdrs.csv", row.names = F)



