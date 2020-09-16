#### THE CODE BELOW IS FROM KATE R'S GITHUB REPO: kricherson/sst_matching.r
library(tidyverse)
library(raster)
library(lubridate)
library(ncdf4)

#Download some example data: Yelloweye rockfish encountered by the NWFSC trawl survey in 2016-2017
yeye_raw<-read_csv("https://www.nwfsc.noaa.gov/data/api/v1/source/trawl.catch_fact/selection.csv?filters=field_identified_taxonomy_dim$scientific_name=Sebastes%20ruberrimus,date_dim$year%3E=2016,date_dim$year%3C=2017&variables=date_yyyymmdd,field_identified_taxonomy_dim$scientific_name,longitude_dd,latitude_dd") 

#What years do we need data for?
years <- unique(year(ymd(yeye_raw$date_yyyymmdd)))

#We're getting SST data from the NOAA OI SST V2 High Resolution Dataset at https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html

#SST file names to download
files <- paste0("sst.day.mean.",years,".nc")

#Download SST data for selected years. This will take a while! # js macbook ~13min
lapply(files, function(filenames){
  download.file(paste0("ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/",filenames), destfile = filenames)
})

#Optional: bounds for cropping the SST data. This makes extracting faster, but note that cropping itself is pretty slow. 
crop_bounds <- extent(min(yeye_raw$longitude_dd %% 360)-0.5, 
                      max(yeye_raw$longitude_dd %% 360)+0.5, 
                      min(yeye_raw$latitude_dd)-0.5, 
                      max(yeye_raw$latitude_dd)+0.5) #This is (xmin, xmax, ymin, ymax)

#Put all of the SST data into a raster stack and crop it. ~2min on JS macbook
sst_stack <- lapply(files,stack) %>% 
  stack(.) %>% 
  crop(., crop_bounds)

#Do a little data manipulation, then match SST to YEYE observations
yeye_dat <- yeye_raw %>% 
  rename(species = `field_identified_taxonomy_dim$scientific_name`)  %>% 
  mutate(longitude = longitude_dd %% 360,
         date = ymd(date_yyyymmdd),
         sst_date = paste0("X", paste(year(date), format.Date(date, "%m"), format.Date(date, "%d"), sep=".")), # this just makes a column that has cell values corresponding to the names of sst_stack columns
         latlon = map2(longitude, latitude_dd, ~ cbind(.x, .y)),
         sst = map2(latlon, sst_date, ~ raster::extract(subset(sst_stack,.y), .x ,method="bilinear")),
         sst = unlist(sst))

#Plot for fun
yeye_sst_plot <- ggplot() + 
  geom_polygon(data = map_data("usa"), aes(x=long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-125.7, -121.0),  ylim = c(36, 49), ratio = 1.3) +
  geom_point(data = yeye_dat, aes(x = longitude_dd, y = latitude_dd, color= sst))