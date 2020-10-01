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
dat$depth <- scale(log(dat$depth))

# load best model
m_po2 <- readRDS("output/wc/model_13_MI.rds") # or use model 11 with o2 rather than po2
# update model with full formula, including breakpoint syntax
#m_po2$formula <- cpue_kg_km2 ~ depth + I(depth^2) + as.factor(year) + temp + breakpt(po2) - 1

# create new data with increments over the range of covariate values observed
nd_temp <- data.frame(temp = seq(min(dat$temp), max(dat$temp), length.out = 100), 
                      depth = 0,
                      po2 = 0,
                      year = 2010L)
nd_po2 <- data.frame(po2 = seq(min(dat$po2), max(dat$po2), length.out = 300), 
                     depth = 0,
                     temp = 0,
                     year = 2010L)
nd_depth <- data.frame(depth = seq(min(dat$depth), max(dat$depth), length.out = 100), 
                      temp = 0,
                      po2 = 0,
                      year = 2010L)

# predict to new data
p_temp <- predict(m_po2, newdata = nd_temp, se_fit = TRUE, re_form = NA, xy_cols = c("longitude", "latitude"))
p_o2 <- predict(m_po2, newdata = nd_po2, se_fit = TRUE, re_form = NA, xy_cols = c("longitude", "latitude"))
p_depth <- predict(m_po2, newdata = nd_depth, se_fit = TRUE, re_form = NA, xy_cols = c("longitude", "latitude"))

# plot predictions with uncertainty
z <- 1.645 # for 90% CI
plot_temp <- ggplot(p_temp, aes(back.convert(temp, attr(dat$temp, "scaled:center"), attr(dat$temp, "scaled:scale")), exp(est), 
                   ymin = exp(est - z * est_se), ymax = exp(est + z * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4) + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks=seq(3,15,3)) +
  labs(x = "Temperature (°C)", y = NULL)
plot_o2 <- ggplot(p_o2, aes(back.convert(po2, attr(dat$po2, "scaled:center"), attr(dat$po2, "scaled:scale")), exp(est), 
                 ymin = exp(est - z * est_se), ymax = exp(est + z * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4) + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Partial Pressure of Oxygen", y = bquote('Population Density'~(kg~km^-2)))
plot_depth <- ggplot(p_depth, aes(exp(back.convert(depth, attr(dat$depth, "scaled:center"), attr(dat$depth, "scaled:scale"))), 
                    exp(est), ymin = exp(est - z * est_se), ymax = exp(est + z * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4) + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Bottom Depth (m)", y = NULL)

multipanel <- gridExtra::grid.arrange(plot_o2, plot_temp, plot_depth, nrow = 1, ncol = 3)
ggsave("output/wc/marginal_effects.pdf", multipanel, width = 9, height = 3, units = "in")