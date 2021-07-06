# Install latest version of sdmTMB ----------------------------------------
#devtools::install_github("pbs-assess/sdmTMB")
# This version uses J-Scope estimates of temp and oxygen

# Initialize with packages and functions ------------------------------------------------
rm(list = ls())
source("code/mi_functions.R")
library(sdmTMB)
library(raster)
library(dplyr)
library(sp)
library(gsw)
library(rgdal)
library(sf)
library(ggplot2)


# 1. Load original data
dat <- load_data(spc = "sablefish", constrain_latitude = T, fit.model= T)
# remove tows w/out DO data here
dat <- dplyr::filter(dat, !is.na(o2))


# save this for back-converting
mean_o2 <- mean(dat$o2)
sd_o2 <- sd(dat$o2)


# transform
dat$log_depth_scaled <- scale(log(dat$depth))
dat$log_depth_scaled2 <- scale(log(dat$depth) ^ 2)
dat$jday_scaled <- scale(dat$julian_day)
dat$jday_scaled2 <- scale((dat$julian_day) ^ 2)
dat$X <- dat$longitude
dat$Y <- dat$latitude

dat$temp <- as.numeric(scale(dat$temp))
dat$o2 <- as.numeric(scale(dat$o2))
dat$po2 <- as.numeric(scale(dat$po2))
dat$mi <- as.numeric(scale(dat$mi))

c_spde <-make_mesh(data = dat, xy_cols = c("X", "Y"), n_knots = 250) # choose # knots

#. fit spatio-temporal model

o2_model <-  sdmTMB(formula = o2 ~ -1 + log_depth_scaled + log_depth_scaled2 
                    +  as.factor(year) + jday_scaled + jday_scaled2,
                    data = dat,
                    time = "year", spde = c_spde, anisotropy = TRUE,
                    silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                    control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

#3. Get prediction
o2_pred <- predict(o2_model,
        newdata = dat,
        return_tmb_object = F)


# add prediction onto dataframe
dat$o2_est <- o2_pred$est* sd_o2 + mean_o2
plot(dat$o2, o2_pred$est, type = "p", pch = 21, xlab = "measured", ylab = "fitted")


# 4. now load j-scope

jscope_dat <- load_data_jscope("sablefish", 2010:2015)

# combine the two

combined_dat <- left_join(x = jscope_dat, 
                          y = dat, 
                           by = "trawl_id")

#  remove tows where haul DO has NA
combined_dat <- dplyr::filter(combined_dat, !is.na(o2.y))

# simple (non spatial plot) of jscope do vs. smoothed do

with(combined_dat,plot(o2.x, o2.y, type = "p",
                       pch = 21,
                       bg = "black",
                       xlab = "J-SCOPE",
                       ylab = "Fitted model to trawl O2 data"
)
)


# make a plot of the difference by space and year

combined_dat <- dplyr::rename(combined_dat,longitude = longitude.x, latitude = latitude.x, 
                              o2_jscope = o2.x,
                              o2_trawl = o2.y,
                              year = year.x)

# put this back in "real" units, mg / l
combined_dat$o2_trawl <- combined_dat$o2_trawl * sd_o2 + mean_o2

combined_dat$delta_o2 <- with(combined_dat, o2_jscope - o2_est)
combined_dat$jscope_minus_trawl <- with(combined_dat, o2_jscope - o2_trawl)

# try again to get latitude and longitude

library(rgdal)

# check zone depending on area of application: WC is primarily 10, BC probably 9, Bering 2 or perhaps 3 near shore
dat_spatial <- SpatialPoints(combined_dat[, c("longitude", "latitude")], proj4string=CRS("+proj=utm +zone=10 +datum=WGS84 +units=km")) 
dat_ll <- as.data.frame(spTransform(dat_spatial, CRS("+proj=longlat")))

combined_dat$Lat <- dat_ll$latitude
combined_dat$Lon <- dat_ll$longitude
require(rnaturalearth)
require(rnaturalearthdata) 
require(rnaturalearthhires)
world = ne_countries(scale = "large", returnclass = "sf")
Canada = subset(world, name == "Canada") # if you wanted to include a boundary between the US and Canada
ggplot(data=world) +
  geom_tile(data=combined_dat, aes(x=Lon, y=Lat, col=jscope_minus_trawl)) + 
  facet_wrap(~year) +
  scale_colour_gradient2() +
  geom_sf(size=0.2) + 
#  geom_sf(data=Canada, size=0.1, fill="gray97") +
  coord_sf(xlim=c(-126, -123.0), ylim=c(43, 49)) + # set for Bering Sea
  labs(x="Longitude", y="Latitude") +
  scale_x_continuous(breaks = round(seq(-126, -123, by = 2),1))
