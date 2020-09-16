# figuring out tidync
# https://github.com/ropensci/tidync
# https://github.com/ropensci/tidync/blob/master/vignettes/netcdf-with-tidync.Rmd

library(tidync)

ice_file <- system.file("extdata", "ifremer", "20171002.nc", package = "tidync", mustWork = TRUE)

tidync(ice_file)

tidync(ice_file) %>% activate(quality_flag)

concentration <- tidync(ice_file) 

hf <- concentration %>% 
  hyper_filter(ni = index > 150, 
               nj = dplyr::between(index, 30, 100))
               
hf %>% 
  hyper_tibble() %>% 
  dplyr::filter(!is.na(concentration)) %>% dplyr::distinct(concentration, quality_flag)

hf %>% active()

hyper_transforms(hf)

## jameal tried to implement tidync and failed. for now.


# example
# file <- system.file("extdata", "oceandata", "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc", package = "tidync")
# library(tidync)
# tidync(file) 
# browseVignettes(package = "tidync")

wc_temp_bottom_file <-  here::here(
  "jacox_ROMS_output",
  "temp_bottom_era5_monthly_1980_2010.nc"
)
tidync(wc_temp_bottom_file)
tidync(wc_temp_bottom_file) %>% hyper_dims()
tidync(wc_temp_bottom_file) %>% hyper_grids()
tidync(wc_temp_bottom_file) %>% activate(temp)

str(tidync(wc_temp_bottom_file) %>% hyper_filter())

lat <- ncvar_get( tidync(wc_temp_bottom_file) ,"lat")

tidync(wc_temp_bottom_file) %>%
  hyper_filter(latitude = between(latitude, -78, -75.8), 
               longitude = between(longitude, 165, 171))


(subs <- tidync(wc_temp_bottom_file) %>% hyper_filter(longitude = longitude < 126, latitude = latitude > 46, time = time > 300))
subs %>% hyper_tibble()
