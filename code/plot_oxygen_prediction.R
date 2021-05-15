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
fit_new_po2_model <- T # do you want to re-fit the spatio-temporal model of oxgyen?


# load data
if(!use_jscope) dat <- load_data(spc = "sablefish", constrain_latitude, fit.model)
if(use_jscope) dat <- load_data_jscope(spc = "sablefish", years = years)

# rescale variables
mean.depth <- mean(log(dat$depth))
sd.depth <- sd(log(dat$depth))

dat$log_depth_scaled <- as.numeric(scale(log(dat$depth)))
dat$log_depth_scaled2 <- dat$log_depth_scaled^2
dat$jday_scaled <- as.numeric(scale(dat$julian_day))
dat$jday_scaled2 <- as.numeric(scale(log(dat$julian_day) ^ 2))

mean.po2 <- mean(dat$po2)
std.po2 <- sd(dat$po2)
mean.do <- mean(dat$o2)
std.do <- sd(dat$o2)
mean.mi <- mean(dat$mi)
std.mi <- sd(dat$mi)
dat$X <- dat$longitude
dat$Y <- dat$latitude

# use this to fit a spatio-temporal model of pO2, so must scale it first
dat$po2 <- as.numeric(scale(dat$po2))
if (fit_new_po2_model) {

# Load best modelfor sablefish
m <- readRDS("output/wc/model_8_sablefish.rds")
c_spde <-make_mesh(data = dat, xy_cols = c("X", "Y"), n_knots = 250) # choose # knots


# removed jday_scaled2 because may not have converged
po2_model <-  sdmTMB(formula = po2 ~ -1 + log_depth_scaled + log_depth_scaled2 
                      + jday_scaled + jday_scaled2,
                     data = dat,
                     time = "year", spde = c_spde, anisotropy = TRUE,
                     silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                     control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

saveRDS(po2_model,"output/wc/po2_model.rds")
}

if(!fit_new_po2_model) po2_model <- readRDS("output/wc/po2_model.rds")

#####
# Now work with the wc_grid, create from scratch

if(fit_new_po2_model) {

grid_cells = readxl::read_excel("data/data/Selection Set 2018 with Cell Corners.xlsx")
if(use_jscope) grid_cells <- dplyr::filter(grid_cells, Cent.Lat >=43)

grid_cells = dplyr::mutate(grid_cells,
                           depth_min = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[1]),
                           depth_max = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[2]))

# convert grid_cells to sp object
grid = SpatialPoints(cbind(grid_cells$Cent.Long,grid_cells$Cent.Lat),
                     proj4string = CRS("+proj=longlat +datum=WGS84"))
r = raster::rasterize(x=grid, y = raster(nrow=length(unique(grid_cells$Cent.Lat)),
                                         ncol=length(unique(grid_cells$Cent.Long))))
rasterToPoints(r)

raster = aggregate(r, fact = 2)
raster = projectRaster(raster, crs = "+proj=tmerc +lat_0=31.96 +lon_0=-121.6 +k=1 +x_0=390000 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# create matrix of point data with coordinates and depth from raster
grid = as.data.frame(rasterToPoints(raster))

newproj = paste("+proj=utm +zone=10 ellps=WGS84")
coordinates(grid_cells) <- c("Cent.Long", "Cent.Lat")
proj4string(grid_cells) <- CRS("+proj=longlat +datum=WGS84")
grid_cells <- spTransform(grid_cells, CRS(newproj))
predict_raster = raster(grid_cells, resolution = c(2778,3704), vals = NULL)
## load custom bathymetry raster
bathy_hiRes <- raster("data/data/bathy_clipped")
bathy_hiRes <- bathy_hiRes / 10 # units were originally decimeters, so convert to meters
# aggregate and project bathymetry to survey grid cells, the absolute minimum resolution of the prediction grid
bathy_raster <- projectRaster(bathy_hiRes, predict_raster, crs = newproj, method="bilinear")
# load Cowcod Conservation Areas, not included in trawl survey, and reproject
CCA = rgdal::readOGR('data/data/kv299cy7357.shp')
CCA = sp::spTransform(CCA, sp::CRS(newproj))
# mask CCA from bathymetry raster used for prediction
bathy_raster = raster::mask(bathy_raster, CCA, inverse = TRUE)
# create matrix of point data with coordinates and depth from raster
wc_grid <- as.data.frame(rasterToPoints(bathy_raster))
colnames(wc_grid) = c("X", "Y", "depth")

# scale covariates
wc_grid$log_depth_scaled <- (log(wc_grid$depth * -1) - mean.depth) / sd.depth
wc_grid$log_depth_scaled2 <- wc_grid$log_depth_scaled ^ 2
wc_grid$X <- wc_grid$X/1000
wc_grid$Y <- wc_grid$Y/1000



wc_grid$loc = seq(1,nrow(wc_grid))

# Create dataframe for fitting
df = expand.grid(loc=unique(wc_grid$loc),
                 year = unique(dat$year))
df = left_join(df,wc_grid)
df$jday_scaled = 0
df$jday_scaled2 = 0
saveRDS(df, file = "output/wc_grid_df.rds")
}

if(!fit_new_po2_model) df <- readRDS("output/wc_grid_df.rds")


pred_po2 <- predict(po2_model,
                    newdata = df,
                    return_tmb_object = F)
# convert estimate (which is scaled) to original po2 units
pred_po2$po2 <- pred_po2$est *std.po2 + mean.po2


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

plotfilename<-"plots/po2_map.pdf"
pdf(file = plotfilename, 
    height = 4,
    width = 8
)
ggplot(us_coast_proj) + geom_sf() +
  geom_raster(data = pred_po2, aes(x = X * 1000, y = Y * 1000, fill = po2)) +
  facet_wrap(~year, ncol = 6) +
  scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
  ylim(ylimits[1], ylimits[2]) +
  scale_fill_viridis_c(limits = c(0.0, 0.1), oob = scales::squish) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14))
dev.off()
system2("open", args = c("-a Skim.app", plotfilename))



