rm(list = ls())
#devtools::install_github("pbs-assess/sdmTMB")
source("code/mi_functions.R")
library(sdmTMB)
library(raster)
library(dplyr)
library(sp)
library(gsw)
library(rgdal)
library(future)
library(sf)
require(rnaturalearth)
require(rnaturalearthdata) 
require(rnaturalearthhires)
require(ggplot2)
library(viridis)
library(tweedie)
library(gridExtra)
library(egg)



# Need to load the data to get scaled stuff
spc <- "sablefish"
fit.model <- F # do you want to fit SDM to environmental variables and impute missing values?
years <- 2010:2015 # designate years to use, must be 2010:2015 if using trawl-based data
constrain_latitude <- F # do you want to constraint trawl data N of 43 degrees latitude
compare_sources <- F # do you want to run models comparing trawl and J-SCOPE covariates?
no_depth <- FALSE # Do you want to run models w/out a depth effect?
use_cv = FALSE # specify whether to do cross validation or not
use_AIC = TRUE # specify whether to use AIC
use_jscope <- F # specify whether to only use J-SCOPE based estimates.  Overrides compare_sources and fit.model
sizeclass <- "p2_p3"

# load and scale data -----------------------------------------------------
dat <- load_data(fit.model= F, spc, constrain_latitude = F)

# rescale variables
mean.depth <- mean(dat$depth)
sd.depth <- sd(dat$depth)
dat$log_depth_scaled <- (scale(log(dat$depth)))
dat$log_depth_scaled2 <- dat$log_depth_scaled^2
dat$jday_scaled <- (scale(dat$julian_day))
dat$jday_scaled2 <- (scale(log(dat$julian_day) ^ 2))
dat$temp <- scale(dat$temp)
dat$po2 <- (scale(dat$po2))
dat$X <- dat$longitude
dat$Y <- dat$latitude
dat$year <- as.factor(dat$year)
dat$mi <- scale(dat$mi)
dat$year <- as.factor(dat$year)
dat$cpue_kg_km2 <- dat$cpue_kg_km2 * (dat$p2+dat$p3)

map_data <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")

# crop if you want; not needed:
us_coast <- st_crop(map_data,
                    c(xmin = -126, ymin = 31, xmax = -110, ymax = 50))
us_coast_proj <- sf::st_transform(us_coast, crs = 32610)



# * 1000 b/c we are in UTM km for model fitting:
if(!use_jscope) {
  xlimits = c(282853, 1025581)
  ylimits = c(3549000, 5366000)
}
if(use_jscope) {
  xlimits = c(283853, 459201)
  ylimits = c(4762418, 5366000)
}

best_model <- readRDS("output/wc/model_7_sablefishp2_p3.rds")

fit_sablefish <- predict(best_model,
                                 newdata = dat,
                                 return_tmb_object = F)

# extract tweedie parameters
phi_index <- which(names(best_model$model$par) == "ln_phi")
theta_index <- which(names(best_model$model$par) == "thetaf")

phi <- exp(best_model$model$par[phi_index])
theta <- best_model$model$par[theta_index]
power <- 1 + exp(theta) / (1 + exp(theta))

# get quantiles of residuals
observed_quantiles <- ptweedie(dat$cpue_kg_km2, mu = exp(fit_sablefish$est), phi = phi, power = power)
fit_sablefish$quantile_residual <- qnorm(observed_quantiles, mean = 0, sd = 1)



ggplot(us_coast_proj) + geom_sf() + 
  geom_point(data = fit_sablefish, aes(x = X * 1000, y = Y * 1000, col = quantile_residual), alpha = 0.6, size = 0.5) +
  facet_wrap(~year, ncol = 3) +
  scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
  ylim(ylimits[1], ylimits[2]) + 
  scale_colour_distiller(type = "div", palette = "RdBu", limits = c(-2.0, 2.00),oob = scales::squish,name = "Quantile Residuals") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text = element_text(size = 12))


ggsave(filename = "plots/residual_catch.png", height = 9, width = 6.5, units = "in")

