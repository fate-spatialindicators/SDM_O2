library(ggplot2)
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)

# load handy functions
source("code/mi_functions.R")

spc <- species.2.plot <- "sablefish"

######
#### Get the data all set up as needed:

                    
dat <- load_data(fit.model, spc, constrain_latitude)

# scale variables - used later to unscale!
dat$temp <- scale(dat$temp)
dat$po2 <- scale(dat$po2)
dat$log_depth_scaled <- scale(log(dat$depth))
dat$log_depth_scaled2 <- dat$log_depth_scaled ^2
dat$mi <- scale(dat$mi)

# load best model
m_po2 <- readRDS("output/wc/model_8_sablefish.rds") #




nd_temp <- data.frame(temp = seq(min(dat$temp), max(dat$temp), length.out = 100),
                      log_depth_scaled = 0,
                      log_depth_scaled2 = 0,
                      po2 = 0,
                      year = 2010L)
nd_temp <- convert_class(nd_temp)


nd_po2 <- data.frame(po2 = seq(min(dat$po2), max(dat$po2), length.out = 300), 
                     temp = 0,
                     log_depth_scaled = 0,
                     log_depth_scaled2 = 0,
                     year = 2010L)
nd_po2 <- convert_class(nd_po2)

nd_depth <- data.frame(log_depth_scaled = seq(min(dat$log_depth_scaled), max(dat$log_depth_scaled), length.out = 100), 
                      temp = 0,
                      po2 = 0,
                      year = 2010L)
nd_depth$log_depth_scaled2 <- nd_depth$log_depth_scaled ^2


nd_depth <- convert_class(nd_depth)

# predict to new data
p_temp <- predict(m_po2, newdata = nd_temp, se_fit = TRUE, re_form = NA)
p_o2 <- predict(m_po2, newdata = nd_po2, se_fit = TRUE, re_form = NA)
p_depth <- predict(m_po2, newdata = nd_depth, se_fit = TRUE, re_form = NA)

# plot predictions with uncertainty
z <- 1.645 # for 90% CI
plot_temp <- ggplot(p_temp, aes(back.convert(temp, attr(dat$temp, "scaled:center"), attr(dat$temp, "scaled:scale")), exp(est), 
                                   ymin = exp(est - z * est_se), ymax = exp(est + z * est_se))) +
                      geom_line() + geom_ribbon(alpha = 0.4) + 
  labs(x = "Temperature (°C)", y = NULL) +
  lims(y = c(0, 1500)) +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) 

plot_o2 <- ggplot(p_o2, aes(back.convert(po2, attr(dat$po2, "scaled:center"), attr(dat$po2, "scaled:scale")), exp(est), 
                 ymin = exp(est - z * est_se), ymax = exp(est + z * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4) + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Partial Pressure of Oxygen", y = bquote('Population Density'~(kg~km^-2))) + 
  lims(y = c(0, 1500)) +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) 

plot_depth <- ggplot(p_depth, aes(exp(back.convert(log_depth_scaled, attr(dat$log_depth_scaled, "scaled:center"), attr(dat$log_depth_scaled, "scaled:scale"))), exp(est), 
                    ymin = exp(est - z * est_se), ymax = exp(est + z * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4) + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Bottom Depth (m)", y = NULL) +
  lims(y = c(0, 1500)) +
theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) 

 
multipanel <- gridExtra::grid.arrange(plot_o2, plot_temp, plot_depth, nrow = 1, ncol = 3)
ggsave("plots/marginal_effects.png", multipanel, width = 8, height = 3, units = "in", device = "png")
