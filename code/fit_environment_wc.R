## Code to fit spatio-temporal model of environmental characteristics 

# Install latest version of sdmTMB ----------------------------------------
devtools::install_github("pbs-assess/sdmTMB")

# Initialize with packages and functions ------------------------------------------------
rm(list = ls())
source("code/mi_functions.R")
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)
library(ggplot2)
library(viridis)

# load and scale data -----------------------------------------------------


dat <- load_data()
spde <- make_spde(x = dat$longitude, y = dat$latitude, n_knots = 250) # choose # knots
mi_model <- sdmTMB(
  formula = mi ~ 0 + as.factor(year),
  data = dat,
  time = "year",
  reml = TRUE,
  spde = spde,
  family = tweedie(link = "log"),
  anisotropy = TRUE,
  spatial_only = F
)

ggplot(predict(mi_model), aes(longitude, latitude, col = exp(est))) + scale_colour_gradientn(colours = plasma(10), limits = c(0, 5)) +
  geom_point(alpha=0.6, size = 0.5) + facet_wrap(~year) + labs(col = "Fitted MI")


o2_model <- try(sdmTMB(
  formula = o2 ~ 0 + as.factor(year),
  data = dat,
  time = "year",
  reml = TRUE,
  spde = spde,
  family = gaussian(link = "log"),
  anisotropy = TRUE,
  spatial_only = FALSE
))

ggplot(predict(o2_model), aes(longitude, latitude, col = est)) + scale_colour_gradientn(colours = plasma(10)) +
  geom_point(alpha=0.6, size = 0.5) + facet_wrap(~year) + labs(col = "Fitted o2")


po2_model <- try(sdmTMB(
  formula = o2 ~ 0 + as.factor(year),
  data = dat,
  time = "year",
  reml = TRUE,
  spde = spde,
  family = gaussian(link = "log"),
  anisotropy = TRUE,
  spatial_only = FALSE
))

ggplot(predict(po2_model), aes(longitude, latitude, col = est)) + scale_colour_gradientn(colours = plasma(10)) +
  geom_point(alpha=0.6, size = 0.5) + facet_wrap(~year) + labs(col = "Fitted po2")

# save to output folder to save time
saveRDS(file = "output/mi_model.rds", mi_model)
saveRDS(file = "output/o2_model.rds", o2_model)
saveRDS(file = "output/po2_model.rds", po2_model)