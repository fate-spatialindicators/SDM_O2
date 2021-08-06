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
fit_new_po2_model <- F # do you want to re-fit the spatio-temporal model of oxgyen?
years.2.plot <- c(2010:2015)

# load data
# trawl insitu-based data

dat <- load_data(spc = "sablefish", constrain_latitude, fit.model)
dat <- dplyr::filter(dat, year == 2010)

nemuro_dat <- load_data_nemuro(spc = "sablefish", constrain_latitude, fit.model)
nemuro_dat <- dplyr::filter(nemuro_dat, year == 2010)

alldat <- left_join(dat, nemuro_dat, by = "trawl_id")

alldat$delta_po2 <- with(alldat, po2.x - po2.y)
alldat$delta_temp <- with(alldat, temp.x - temp.y)
alldat$delta_mi <- with(alldat, mi.x - mi.y)

alldat <- dplyr::filter(alldat, !is.na(delta_po2), !is.na(delta_temp), !is.na(delta_mi))





# sean's code
map_data <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")

# crop if you want; not needed:
us_coast <- st_crop(map_data,
                    c(xmin = -126, ymin = 31, xmax = -110, ymax = 50))
us_coast_proj <- sf::st_transform(us_coast, crs = 32610)



# * 1000 b/c we are in UTM km for model fitting:

  xlimits = c(282853, 1025581)
  ylimits = c(3549000, 5366000)


tempmap <-ggplot(us_coast_proj) + geom_sf() +
  geom_point(data = alldat, aes(x = longitude.x * 1000, y = latitude.x * 1000, col = delta_temp)) +
  scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
  ylim(ylimits[1], ylimits[2]) +
  scale_colour_viridis_c(limits = c(-2, 1), oob = scales::squish,name = bquote(Delta ~"Temp")) +
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

temp_compare <- ggplot(alldat) + 
  geom_point(aes(x = temp.x, y = temp.y, col = depth.x)) +
  scale_colour_viridis_c(name = "Depth") +
  labs(x = "In situ Temperature", y = "NEMURO Temperature") +
  theme_bw() +
  theme(panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text = element_text(size = 12)) +
  geom_abline(intercept = 0, slope = 1)
  

grid.arrange(tempmap, temp_compare, ncol = 2)
ggsave("plots/compare_temperature.png")


po2map <-ggplot(us_coast_proj) + geom_sf() +
  geom_point(data = alldat, aes(x = longitude.x * 1000, y = latitude.x * 1000, col = delta_po2)) +
  scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
  ylim(ylimits[1], ylimits[2]) +
  scale_colour_viridis_c(limits = c(-0.1, 0.05), oob = scales::squish,name = bquote(Delta ~pO[2])) +
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

po2_compare <- ggplot(alldat) + 
  geom_point(aes(x = po2.x, y = po2.y, col = depth.x)) +
  scale_colour_viridis_c(name = "Depth") +
  labs(x = "In situ pO2", y = "NEMURO pO2") +
  theme_bw() +
  theme(panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text = element_text(size = 12)) +
  geom_abline(intercept = 0, slope = 1)


grid.arrange(po2map, po2_compare, ncol = 2)
ggsave("plots/compare_po2.png")
mimap <-ggplot(us_coast_proj) + geom_sf() +
  geom_point(data = alldat, aes(x = longitude.x * 1000, y = latitude.x * 1000, col = delta_mi)) +
  scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
  ylim(ylimits[1], ylimits[2]) +
  scale_colour_viridis_c(limits = c(-3, 1.5), oob = scales::squish,name = bquote(Delta ~"MI")) + 
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

mi_compare <- ggplot(alldat) + 
  geom_point(aes(x = mi.x, y = mi.y, col = depth.x)) +
  scale_colour_viridis_c(name = "Depth") +
  labs(x = "In situ MI", y = "NEMURO MI") +
  theme_bw() +
  theme(panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text = element_text(size = 12)) +
  geom_abline(intercept = 0, slope = 1)


grid.arrange(mimap, mi_compare, ncol = 2)
ggsave("plots/compare_mi.png")
