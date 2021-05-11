rm(list = ls())
library(ggplot2)
library(raster)
library(rasterize)
library(sp)

library(sdmTMB)
library(dplyr)

grid_cells = readxl::read_excel("data/data/Selection Set 2018 with Cell Corners.xlsx")
grid_cells = dplyr::mutate(grid_cells,
                           depth_min = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[1]),
                           depth_max = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[2]))
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
newproj = paste("+proj=utm +zone=10 +ellps=WGS84")
coordinates(grid_cells) <- c("Cent.Long", "Cent.Lat")
proj4string(grid_cells) <- CRS("+proj=longlat +ellips=WGS84")
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
wc_grid$X <- wc_grid$X / 1000
wc_grid$Y <- wc_grid$Y / 1000
saveRDS(wc_grid, file=paste0("data/wc_grid.rds")) # save prediction grid

###### 
wc_grid <- readRDS("survey_data/wc_grid.rds")
