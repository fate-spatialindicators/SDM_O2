library(ggplot2)
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)

# load handy functions
source("code/mi_functions.R")

dat <- load_data()

# scale variables
dat$temp <- scale(dat$temp)
dat$po2 <- scale(dat$po2)

# load best model
m_po2 <- readRDS("output/wc/model_13_MI.rds")

# create new data with increments over the range of covariate values observed
nd_temp <- data.frame(temp = seq(min(dat$temp), max(dat$temp), length.out = 100), 
                      depth = 0,
                      po2 = 0,
                      year = 2010L)
nd_po2 <- data.frame(po2 = seq(min(dat$po2), max(dat$po2), length.out = 100), 
                     depth = 0,
                     temp = 0,
                     year = 2010L)

# BELOW 2 LINES CAUSING R TO CRASH
p_temp <- predict(m_po2, newdata = nd_temp, se_fit = TRUE, re_form = NA)
p_o2 <- predict(m_po2, newdata = nd_po2, se_fit = TRUE, re_form = NA)

ggplot(p_temp, aes(back.convert(temp, 6.832134, 2.007563), exp(est), ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4)
ggplot(p_o2, aes(back.convert(po2, 0.03841377, 0.03037818), exp(est), ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4)