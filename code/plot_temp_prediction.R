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
require(rnaturalearth)
require(rnaturalearthdata) 
require(rnaturalearthhires)
require(ggplot2)
library(viridis)
library(sf)
library(gridExtra)
library(egg)


# load handy functions
source("code/mi_functions.R")

# Need to load the data to get scaled stuff
fit.model <- F # do you want to fit SDM to environmental variables and impute missing values?
years <- 2010:2015 # designate years to use, must be 2010:2015 if using trawl-based data
constrain_latitude <- F # do you want to constraint trawl data N of 43 degrees latitude
compare_sources <- F # do you want to run models comparing trawl and J-SCOPE covariates?
no_depth <- FALSE # Do you want to run models w/out a depth effect?
use_cv = FALSE # specify whether to do cross validation or not
use_AIC = TRUE # specify whether to use AIC
use_jscope <- F # specify whether to only use J-SCOPE based estimates.  Overrides compare_sources and fit.model
fit_new_temp_model <- F # do you want to re-fit the spatio-temporal model of oxgyen?
years.2.plot <- c(2010:2015)

# load data
if(!use_jscope) dat <- load_data(spc = "sablefish", constrain_latitude, fit.model)
if(use_jscope) dat <- load_data_jscope(spc = "sablefish", years = years)

# rescale variables
dat$log_depth_scaled <- (scale(log(dat$depth)))
dat$log_depth_scaled2 <- dat$log_depth_scaled^2
dat$jday_scaled <- (scale(dat$julian_day))
dat$jday_scaled2 <- (scale(log(dat$julian_day) ^ 2))
dat$X <- dat$longitude
dat$Y <- dat$latitude
dat$year <- as.factor(dat$year)
dat$temp <- scale(dat$temp)

if (fit_new_temp_model) {

c_spde <-make_mesh(data = dat, xy_cols = c("X", "Y"), n_knots = 250) # choose # knots
temp_model <-  sdmTMB(formula = temp ~ -1 + log_depth_scaled + log_depth_scaled2 
                      + jday_scaled + jday_scaled2 + year,
                     data = dat,
                     time = "year", spde = c_spde, anisotropy = TRUE,
                     silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                     control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

saveRDS(temp_model,"output/wc/temp_model.rds")
}

if(!fit_new_temp_model) temp_model <- readRDS("output/wc/temp_model.rds")

#####
# Now work with the wc_grid
df <- readRDS("output/wc_grid_df.rds")
df$year <- as.factor(df$year)
df$log_depth_scaled <- as(df$log_depth_scaled, Class = "matrix")
df$log_depth_scaled2 <- as(df$log_depth_scaled2, Class = "matrix")
df$jday_scaled <- as(df$jday_scaled, Class = "matrix")
df$jday_scaled2 <- as(df$jday_scaled2, Class = "matrix")
pred_temp <- predict(temp_model,
                    newdata = df,
                    return_tmb_object = F)
# convert estimate (which is scaled) to original temp units
pred_temp$temp <- back.convert(pred_temp$est, attr(dat$temp,"scaled:center"), attr(dat$temp, "scaled:scale"))

pred_temp <- dplyr::filter(pred_temp, year %in% years.2.plot)
if(!use_jscope) saveRDS(pred_temp, file = "output/wc/pred_temp.rds")
if(use_jscope) saveRDS(pred_temp, file = "output/wc/pred_temp_jscope.rds")
# sean's code
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


ggplot(us_coast_proj) + geom_sf() +
  geom_raster(data = pred_temp, aes(x = X * 1000, y = Y * 1000, fill = temp)) +
  facet_wrap(~year, ncol = 3) +
  scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
  ylim(ylimits[1], ylimits[2]) +
  scale_fill_viridis_c(limits = c(3, 8), oob = scales::squish, name = "Temperature (C)") +
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

ggsave(filename = "plots/tempmap.png", height = 9, width = 6.5, units = "in")