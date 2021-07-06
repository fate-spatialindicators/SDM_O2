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

library(gridExtra)
library(egg)



# Need to load the data to get scaled stuff
fit.model <- F # do you want to fit SDM to environmental variables and impute missing values?
years <- 2010:2015 # designate years to use, must be 2010:2015 if using trawl-based data
constrain_latitude <- F # do you want to constraint trawl data N of 43 degrees latitude
compare_sources <- F # do you want to run models comparing trawl and J-SCOPE covariates?
no_depth <- FALSE # Do you want to run models w/out a depth effect?
use_cv = FALSE # specify whether to do cross validation or not
use_AIC = TRUE # specify whether to use AIC
use_jscope <- F # specify whether to only use J-SCOPE based estimates.  Overrides compare_sources and fit.model
years.2.plot <- c(2010:2015)

# load data
if(!use_jscope) dat <- load_data(spc = "sablefish", constrain_latitude, fit.model)
if(use_jscope) dat <- load_data_jscope(spc = "sablefish", years = years)

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

best_model <- readRDS("output/wc/model_8_sablefish.rds")
pred_po2 <- readRDS("output/wc/pred_po2.rds")
pred_temp <-readRDS("output/wc/pred_temp.rds")

df <- readRDS("output/wc_grid_df.rds")
df$year <- as.factor(df$year)
df$log_depth_scaled <- as(df$log_depth_scaled, Class = "matrix")
df$log_depth_scaled2 <- as(df$log_depth_scaled2, Class = "matrix")
df$jday_scaled <- as(df$jday_scaled, Class = "matrix")
df$jday_scaled2 <- as(df$jday_scaled2, Class = "matrix")
df$po2 <- as(pred_po2$est, Class = "matrix")
df$temp <- as(pred_temp$est, Class = "matrix")


predict_sablefish <- predict(best_model,
                                 newdata = df,
                                 return_tmb_object = F)

df_highpo2 <- df
po2_threshold <- best_model$model$par[length(best_model$model$par)]
df_highpo2$po2[df_highpo2$po2<=po2_threshold] = po2_threshold

predict_sablefish_highpo2 <-predict(best_model,
                                    newdata = df_highpo2,
                                    return_tmb_object = F)
delta_predict <- predict_sablefish
delta_predict$delta <- ((predict_sablefish$est)) - ((predict_sablefish_highpo2$est))

ggplot(us_coast_proj) + geom_sf() +
  geom_raster(data = delta_predict, aes(x = X * 1000, y = Y * 1000, fill = delta)) +
  facet_wrap(~year, ncol = 3) +
  scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
  ylim(ylimits[1], ylimits[2]) +
  scale_fill_viridis_c(limits = c(-1.5, 0), oob = scales::squish, name = "log Effect Size") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text = element_text(size = 12))

ggsave("plots/spatial_effects.png", width = 6.5, height = 9, units = "in", device = "png")

saveRDS(delta_predict, file = "output/wc/delta_predict.rds")
