rm(list = ls())
devtools::install_github("pbs-assess/sdmTMB")
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
fit_new_po2_model <- T # do you want to re-fit the spatio-temporal model of oxgyen?
years.2.plot <- c(2010:2015)

# load data
dat <- load_data(spc = "sablefish", constrain_latitude, fit.model)


# rescale variables
dat$X <- dat$longitude
dat$Y <- dat$latitude
dat$year <- as.factor(dat$year)



# sean's code
map_data <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")

# crop if you want; not needed:
us_coast <- st_crop(map_data,
                    c(xmin = -126, ymin = 31, xmax = -110, ymax = 50))
us_coast_proj <- sf::st_transform(us_coast, crs = 32610)



# * 1000 b/c we are in UTM km for model fitting:

  xlimits = c(282853, 1025581)
  ylimits = c(3549000, 5366000)

red.dat <- dplyr::filter(dat, year %in% c(2011, 2015))

po2_plot <- ggplot(us_coast_proj) + geom_sf() +
  geom_point(data = red.dat, aes(x = X * 1000, y = Y * 1000, col = po2), size = 0.1, alpha = 1.0) +
  facet_wrap(~year, ncol = 1) +
  scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
  ylim(ylimits[1], ylimits[2]) +
  scale_colour_viridis_c(limits = c(0, 0.1), oob = scales::squish, name = bquote(pO[2]), breaks = c(0, 0.05, 0.1)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        , strip.background = element_blank()
        , strip.text = element_blank()
  ) + 
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.position = "bottom") + 
  guides(colour = guide_colourbar(title.position="top", title.hjust=0.5))


temp_plot <- ggplot(us_coast_proj) + geom_sf() +
  geom_point(data = red.dat, aes(x = X * 1000, y = Y * 1000, col = temp), size = 0.1, alpha = 1.0) +
  facet_wrap(~year, ncol = 1) +
  scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
  ylim(ylimits[1], ylimits[2]) +
  scale_colour_viridis_c(limits = c(3, 11), oob = scales::squish, name = "Temperature (C)") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        , strip.background = element_blank()
        , strip.text = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) +
  theme(axis.title.y= element_text(colour = "white")) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.position = "bottom") + 
  guides(colour = guide_colourbar(title.position="top", title.hjust=0.5))
         

mi_plot <-  ggplot(us_coast_proj) + geom_sf() +
  geom_point(data = red.dat, aes(x = X * 1000, y = Y * 1000, col = mi), size = 0.1, alpha = 1.0) +
  facet_wrap(~year, ncol = 1) +
  scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
  ylim(ylimits[1], ylimits[2]) +
  scale_colour_viridis_c(limits = c(0, 6), oob = scales::squish, name = "Metabolic Index") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        , strip.background = element_blank()
        , strip.text = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) +
  theme(axis.title.y= element_text(colour = "white")) +
  theme(legend.text = element_text(size = 12)) + 
  theme(legend.position = "bottom") + 
  guides(colour = guide_colourbar(title.position="top", title.hjust=0.5))

p <- grid.arrange(po2_plot, temp_plot, mi_plot, ncol = 3)

ggsave(filename = "plots/env_summary.pdf",device = "pdf", p, height = 6, width = 9, units = "in")

# this creates a lot of white space in between columns, fixing this in illustrator
