#devtools::install_github("pbs-assess/sdmTMB")

source("code/mi_functions.R")
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)
library(rgdal)
library(future)
library(ggplot2)

library(viridis)
cols <- viridis(2)

dat_bc <- load_bc_data(species.name = "sablefish", survey.name = "SYN WCHG")
dat_wc <- load_data()
dat_bc$year <- as.factor(dat_bc$year)
dat_bc$region <- "BC"
dat_wc$region <- "WC"
dat_wc <- dat_wc %>%
  select(species, year, cpue_kg_km2, o2, temp, depth, mi, po2, longitude, latitude, region)

dat <- rbind(dat_wc, dat_bc)


hist(dat_bc$po2, main = "BC")
hist(dat_wc$po2, main = "WC")
ggplot(dat, aes(x = depth, y = po2,col = region)) + geom_point() +
  scale_color_manual(values=cols)
cols = plasma(6)
ggplot(dat_bc, aes(x = depth, y = po2, col = year)) + geom_point() +
  scale_color_manual(values=cols)
