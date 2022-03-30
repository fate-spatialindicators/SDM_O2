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



# load data
dat <- readRDS("survey_data/joined_nwfsc_data.rds")
dat_by_size <- readRDS("data/survey_data/sablefish_size_dist.rds")
dat = dplyr::filter(dat, species == "sablefish", year%in%seq(2010,2015))
dat <- left_join(dat, dat_by_size, by = "trawl_id")
# remove tows where there was positive catch but no length measurements
dat <- dplyr::filter(dat, !is.na(p1))

# get julian day
dat$julian_day <- rep(NA, nrow(dat))
for (i in 1:nrow(dat)) dat$julian_day[i] <- as.POSIXlt(dat$date[i], format = "%Y-%b-%d")$yday



dat$date <- as.POSIXct(dat$date, format = "%Y-%b-%d", tz = "")
dat$datex <- as.Date(dat$date)
# compute metabolic index (mi) --------------------------------------------
# converted from Halle Berger matlab script

#O2 from trawl data is in ml/l 
# just in case, remove any missing or nonsense values from sensors
dat <- dplyr::filter(dat, !is.na(o2), !is.na(sal), !is.na(temp), is.finite(sal))
dat <- calc_po2_mi(dat)
dat <- dplyr::filter(dat, !is.na(temp), !is.na(mi))

# prepare data and models -------------------------------------------------

dat <- dplyr::select(dat, trawl_id, species, year, longitude_dd, latitude_dd, cpue_kg_km2,
                     o2, temp, depth, mi, po2, julian_day, pass, p1, p2, p3, p4, datex)




shallow.dat <- dplyr::filter(dat,depth < 200)
shallow.dat$year <- as.factor(shallow.dat$year)
multipanel <- ggplot() + 
  geom_point(data = shallow.dat, aes(x = datex, y = po2, col = latitude_dd)) + 
  facet_wrap(~year, ncol = 3, scales = "free") +
  scale_colour_viridis_c(name = "Latitude") +
  scale_x_date(date_labels = "%b %d") +
  coord_cartesian(ylim = c(0, 16.0)) + 
  scale_y_continuous(expand = c(0,0))+
  labs(x = "Date", y = bquote(pO[2]~"(kPa)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        , panel.ontop = element_blank()
        , strip.background= element_blank()
        , strip.placement= "inside"
        , plot.title.position= "plot"
        , strip.text = element_text(size = 12, hjust = 0, vjust = 0)
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text = element_text(size = 12))
ggsave("plots/shelf_po2.png", multipanel, width = 10, height = 5, units = "in", device = "png")

