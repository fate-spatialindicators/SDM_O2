#devtools::install_github("pbs-assess/sdmTMB")
library(ggplot2)
library(raster)
library(rasterize)
library(sp)

library(sdmTMB)
library(dplyr)

# haul data includes environmental covariates with location information
haul = nwfscSurvey::PullHaul.fn(SurveyName = "NWFSC.Combo")

# read in the grid cell data from the survey design
grid_cells = readxl::read_excel("data/data/Selection Set 2018 with Cell Corners.xlsx")
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

# Figure out the grid cell corresponding to each tow location
haul$Cent.Lat = NA
haul$Cent.Lon = NA
haul$Cent.ID = NA
for(i in 1:nrow(haul)) {
  indx = which(grid_cells$NW.LAT > haul$latitude_dd[i] &
      grid_cells$SW.LAT < haul$latitude_dd[i] &
      grid_cells$NW.LON < haul$longitude_dd[i] &
      grid_cells$NE.LON > haul$longitude_dd[i])
  if(length(indx) > 0) {
    haul$Cent.ID[i] = grid_cells$Cent.ID[indx]
    haul$Cent.Lat[i] = grid_cells$Cent.Lat[indx]
    haul$Cent.Lon[i] = grid_cells$Cent.Long[indx]
  }
}

# project lat/lon to UTM, after removing missing values and unsatisfactory hauls
haul = haul %>% filter(!is.na(Cent.Lon), performance == "Satisfactory")

haul_trans = haul
coordinates(haul_trans) <- c("Cent.Lon", "Cent.Lat")
proj4string(haul_trans) <- CRS("+proj=longlat +datum=WGS84")
newproj = paste("+proj=utm +zone=10 ellps=WGS84")
haul_trans <- spTransform(haul_trans, CRS(newproj))
haul_trans = as.data.frame(haul_trans)
haul_trans$Cent.Lon = haul_trans$Cent.Lon/10000
haul_trans$Cent.Lat = haul_trans$Cent.Lat/10000
haul_trans$year = as.numeric(substr(haul_trans$date_yyyymmdd,1,4))

haul$X = haul_trans$Cent.Lon
haul$Y = haul_trans$Cent.Lat
haul$year = haul_trans$year
#haul$year_centered = haul$year - mean(unique(haul$year))

# center and scale depth, removing NAs
haul = dplyr::filter(haul, !is.na(depth_hi_prec_m))
haul$log_depth_scaled = scale(log(haul$depth_hi_prec_m))
haul$log_depth_scaled2 = haul$log_depth_scaled ^ 2

coordinates(grid_cells) <- c("Cent.Long", "Cent.Lat")
proj4string(grid_cells) <- CRS("+proj=longlat +datum=WGS84")
grid_cells <- spTransform(grid_cells, CRS(newproj))

# make prediction raster roughly from grid_cell centroids, given standard cell dimensions (here in meters, converted from nm)
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
wc_grid$log_depth_scaled <- (log(wc_grid$depth * -1) - mean(log(haul$depth_hi_prec_m))) / sd(log(haul$depth_hi_prec_m))
wc_grid$log_depth_scaled2 <- wc_grid$log_depth_scaled ^ 2
wc_grid$X <- wc_grid$X/10000
wc_grid$Y <- wc_grid$Y/10000

saveRDS(wc_grid, file=paste0("data/wc_grid.rds")) # save prediction grid

n_knots = 300

o2haul = dplyr::filter(haul, year>=2009, !is.na(o2_at_gear_ml_per_l_der))

# try to add in seasonal component
o2haul$month = as.numeric(substr(o2haul$date_yyyymmdd, 5,6))
o2haul$day = as.numeric(substr(o2haul$date_yyyymmdd, 7,8))
o2haul$jday = date::mdy.date(o2haul$month, o2haul$day, o2haul$year) -
  date::mdy.date(1, 1, o2haul$year)
o2haul$jday_scaled = scale(o2haul$jday)
o2haul$jday_scaled2 = o2haul$jday_scaled ^ 2

# fit first model to raw geospatial data
c_spde <- make_spde(o2haul$X, o2haul$Y, n_knots = n_knots)
o2_model <- sdmTMB(formula = o2_at_gear_ml_per_l_der ~ -1 + log_depth_scaled +
    log_depth_scaled2 + as.factor(year) + jday_scaled + jday_scaled2,
      data = o2haul,
      time = "year", spde = c_spde, anisotropy = TRUE,
      silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
      control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

# load prediction grid
wc_grid = readRDS(paste0("data/wc_grid.rds"))
wc_grid$loc = seq(1,nrow(wc_grid))

# truncate limits based on haul filters for OR above
df = expand.grid(loc=unique(wc_grid$loc),
  year = unique(o2haul$year))

df = left_join(df,wc_grid)
df$jday_scaled = 0
df$jday_scaled2 = 0
df$log_depth_scaled = matrix(df$log_depth_scaled, ncol=1)
df$log_depth_scaled2 = matrix(df$log_depth_scaled2, ncol=1)
df$jday_scaled = matrix(df$jday_scaled, ncol=1)
df$jday_scaled2 = matrix(df$jday_scaled2, ncol=1)

pred_o2 = predict(o2_model,
  newdata=df,
  xy_cols=c("X","Y"),
  return_tmb_object = TRUE)

save(o2_model,pred_o2,file="wc_o2.Rdata")


# Do the same model with temperature

# try to add in seasonal component
haul$month = as.numeric(substr(haul$date_yyyymmdd, 5,6))
haul$day = as.numeric(substr(haul$date_yyyymmdd, 7,8))
haul$jday = date::mdy.date(haul$month, haul$day, haul$year) -
  date::mdy.date(1, 1, haul$year)
haul$jday_scaled = scale(haul$jday)
haul$jday_scaled2 = haul$jday_scaled ^ 2

# fit first model to raw geospatial data
c_spde <- make_spde(haul$X, haul$Y, n_knots = n_knots)
temp_model <- sdmTMB(formula = temperature_at_gear_c_der ~ -1 + log_depth_scaled +
    log_depth_scaled2 + as.factor(year) + jday_scaled + jday_scaled2,
  data = haul,
  time = "year", spde = c_spde, anisotropy = TRUE,
  silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,#family = tweedie(link = "log"),
  control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

# load prediction grid
wc_grid = readRDS(paste0("data/wc_grid.rds"))
wc_grid$loc = seq(1,nrow(wc_grid))

# truncate limits based on haul filters for OR above
df = expand.grid(loc=unique(wc_grid$loc),
  year = unique(haul$year))
df = left_join(df,wc_grid)
df$jday_scaled = 0
df$jday_scaled2 = 0
df$log_depth_scaled = matrix(df$log_depth_scaled, ncol=1)
df$log_depth_scaled2 = matrix(df$log_depth_scaled2, ncol=1)
df$jday_scaled = matrix(df$jday_scaled, ncol=1)
df$jday_scaled2 = matrix(df$jday_scaled2, ncol=1)

pred_temp = predict(temp_model,
  newdata=df,
  xy_cols=c("X","Y"),
  return_tmb_object = TRUE)

save(temp_model,pred_temp,file="wc_temp.Rdata")



