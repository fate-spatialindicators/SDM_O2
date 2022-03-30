rm(list = ls())
source("code/mi_functions.R")
library(raster)
library(dplyr)
library(sp)
library(gsw)
library(rgdal)
require(rnaturalearth)
require(rnaturalearthdata) 
require(rnaturalearthhires)
require(ggplot2)
library(viridis)
library(sf)
library(gridExtra)
library(egg)



# Set Run Specifications --------------------------------------------------
spc <- "sablefish"
fit.model <- F # do you want to fit SDM to environmental variables and impute missing values?
years <- 2010:2015 # designate years to use, must be 2010:2015 if using trawl-based data
constrain_latitude <- F # do you want to constraint trawl data N of 43 degrees latitude
compare_sources <- F # do you want to run models comparing trawl and J-SCOPE covariates?
no_depth <- FALSE # Do you want to run models w/out a depth effect?
use_jscope <- F # specify whether to only use J-SCOPE based estimates.  Overrides compare_sources and fit.model
sizeclass <- "p2_p3"
k_folds <- 5
constrain_depth <- F  # here you can choose to only look at depths within a certain range. 

# load and scale data -----------------------------------------------------
dat <- load_data(fit.model= F, spc, constrain_latitude = F)
# to test, remove depths > 800 m
if (constrain_depth) dat <- dplyr::filter(dat, depth>=800)

dat$temp <- (scale(dat$temp))
dat$mi <- (scale(dat$mi))
dat$po2 <- (scale(dat$po2))
dat$log_depth_scaled <- scale(log(dat$depth))
dat$log_depth_scaled2 <- with(dat, log_depth_scaled ^ 2)
dat$log_depth_scaled3 <- with(dat, log_depth_scaled^3)
dat$jday_scaled <- scale(dat$julian_day)
dat$jday_scaled2 <- with(dat, jday_scaled ^ 2)
dat$X <- dat$longitude
dat$Y <- dat$latitude
dat$cpue_kg_km2 <- dat$cpue_kg_km2 * (dat$p2+dat$p3)

# make year a factor
dat$year <- as.factor(dat$year)

map_data <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")

# crop if you want; not needed:
us_coast <- st_crop(map_data,
                    c(xmin = -126, ymin = 31, xmax = -110, ymax = 50))
us_coast_proj <- sf::st_transform(us_coast, crs = 32610)


# Plot sablefish catch rates #### 

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
  geom_point(data = dat, aes(x = X * 1000, y = Y * 1000, col = log(cpue_kg_km2)), alpha = 0.6, size = 0.5) +
  facet_wrap(~year, ncol = 3) +
  scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
  ylim(ylimits[1], ylimits[2]) + 
  scale_colour_viridis(limits = c(4, 8),oob = scales::squish,name = bquote('log Catch Rate'~(kg~km^-2))) +
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

ggsave(filename = "plots/sablemap.png", height = 9, width = 6.5, units = "in")