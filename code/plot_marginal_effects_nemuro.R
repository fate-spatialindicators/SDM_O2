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

                    
dat <- load_data_nemuro()

# scale variables - used later to unscale!
dat$temp <- scale(dat$temp)
dat$po2 <- scale(dat$po2)
dat$log_depth_scaled <- scale(log(dat$depth))
dat$log_depth_scaled2 <- dat$log_depth_scaled ^2
dat$mi <- scale(dat$mi)

# load best model
m_mi <- readRDS("output/wc/model_9_sablefish_nemuro.rds") #




nd_mi <- data.frame(mi = seq(min(dat$mi), max(dat$mi), length.out = 100),
                      log_depth_scaled = 0,
                      log_depth_scaled2 = 0,
                      year = 2010L)
nd_mi <- convert_class(nd_mi)


nd_depth <- data.frame(log_depth_scaled = seq(min(dat$log_depth_scaled), max(dat$log_depth_scaled), length.out = 100), 
                      mi = 0,
                      year = 2010L)
nd_depth$log_depth_scaled2 <- nd_depth$log_depth_scaled ^2


nd_depth <- convert_class(nd_depth)

# predict to new data
p_mi <- predict(m_mi, newdata = nd_mi, se_fit = TRUE, re_form = NA)
p_depth <- predict(m_mi, newdata = nd_depth, se_fit = TRUE, re_form = NA)

# plot predictions with uncertainty
z <- 1.645 # for 90% CI
plot_mi <- ggplot(p_mi, aes(back.convert(mi, attr(dat$mi, "scaled:center"), attr(dat$mi, "scaled:scale")), exp(est), 
                                   ymin = exp(est - z * est_se), ymax = exp(est + z * est_se))) +
                      geom_line() + geom_ribbon(alpha = 0.4) + 
  labs(x = "Metabolic Index", y = NULL) +
  scale_y_continuous(limits = c(0, 2250), expand = expansion(mult = c(0, 0.0))) +
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
  scale_y_continuous(limits = c(0, 2250), expand = expansion(mult = c(0, 0.0))) +
  labs(x = "Bottom Depth (m)", y = NULL) +
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

 
multipanel <- gridExtra::grid.arrange(plot_mi, plot_depth, nrow = 1, ncol = 2)
ggsave("plots/marginal_effects_nemuro.png", multipanel, width = 8, height = 3, units = "in", device = "png")
