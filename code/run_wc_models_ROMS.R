# Dependencies: run Extracto-ROMS.R first
devtools::install_github("pbs-assess/sdmTMB")
library(sdmTMB)
library(dplyr)
library(sp)

# get survey data

ROMS.names <- c("ROMS_oxygen_bottom_era5_monthly","ROMS_temp_bottom_era5_monthly")
ROMS.RDS.names <- paste0("survey_data/joined_nwfsc_data",ROMS.names,".rds")

#dat = readRDS("survey_data/joined_nwfsc_data.rds")
dat_temp_bottom <- readRDS(ROMS.RDS.names[2])
dat_oygen_bottom <- readRDS(ROMS.RDS.names[1])
# lapply(ROMS.RDS.names, function(filenames){
#   readRDS(ROMS.RDS.names[i])
# })
dat <- dat_temp_bottom %>%
  left_join(dat_oygen_bottom)

#dat$date <- ymd(dat$date)
#dat$month <- month(dat$date)
# saveRDS(dat,
#         file=here::here(
#           "survey_data",
#           "joined_nwfsc_data_month.rds")
# )



dplyr::group_by(dat, species) %>% 
  summarize(min = min(latitude_dd[which(cpue_kg_km2 > 0)]),
            max = max(latitude_dd[which(cpue_kg_km2 > 0)])) %>% 
  as.data.frame() %>% arrange(min)

# UTM transformation
dat_ll = dat
coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
# convert to utm with spTransform
dat_utm = spTransform(dat_ll, 
                      CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
# convert back from sp object to data frame
dat = as.data.frame(dat_utm)
dat = dplyr::rename(dat, longitude = longitude_dd, 
                    latitude = latitude_dd)

df = expand.grid("species" = unique(dat$species),
                 spatial_only=c(FALSE), 
                 depth_effect = c(TRUE,FALSE),
                 time_varying = c(FALSE),
                 covariate = c("ROMS_temp_bottom_era5_monthly"))#,"ROMS_oxygen_bottom_era5_monthly")

saveRDS(df, "output/wc/jacox_ROMS/models_ROMS.RDS")

start.time <- Sys.time()
for(i in 1:nrow(df)) {
  
  # filter by species, and select range within occurrences
  sub = dplyr::filter(dat, 
                      species == df$species[i])# %>% 
    #dplyr::filter(latitude > min(latitude[which(cpue_kg_km2>0)]),
    #              latitude <= max(latitude[which(cpue_kg_km2>0)]),
    #              longitude > min(longitude[which(cpue_kg_km2>0)]),
    #              longitude < max(longitude[which(cpue_kg_km2>0)]))

  # drop points with missing values and recale to be N(0,1)
  sub = dplyr::filter(sub, 
    !is.na(ROMS_oxygen_bottom_era5_monthly),
    !is.na(ROMS_temp_bottom_era5_monthly),
    !is.na(depth)) %>% 
    dplyr::mutate(depth = scale(log(depth)),
      o2 = scale(ROMS_oxygen_bottom_era5_monthly),
      temp = scale(ROMS_temp_bottom_era5_monthly))
  
  # filter years based on covariate. for ROMS right now, only do 2003-2010
  sub = dplyr::filter(sub, year%in%seq(2003,2010))
  
  # if(df$covariate[i]=="o2") {
  #   sub = dplyr::filter(sub, year%in%seq(2010,2015))
  # } else {
  #   # temp observed for all years
  #   sub = sub
  # }
  
  # rename variables to make code generic
  sub = dplyr::rename(sub, enviro = as.character(df$covariate[i]))
  
  # make spde
  spde <- make_spde(x = sub$longitude, y = sub$latitude, 
                    n_knots = 250) # ew recommends 250 as a minimum, max of 350
  
  formula = paste0("cpue_kg_km2 ~ -1")
  if(df$depth_effect[i]==TRUE) {
    formula = paste0(formula, " + depth + I(depth^2)")
  }
  
  time_formula = "~ -1"
  if(df$time_varying[i]==TRUE) {
    time_formula = paste0(time_formula, " + ", 
                          "enviro", " + I(","enviro","^2)")
    time_varying = as.formula(time_formula)
    time = "year"
  } else {
    formula = paste0(formula, " + ", 
                     "enviro", " + I(","enviro","^2)")
    time_varying = NULL
    time = "year"
  }
  formula = paste0(formula, " + as.factor(year)")
  
  # fit model
  m <- try( sdmTMB( 
    formula = as.formula(formula),
    #formula = cpue_kg_km2 ~ as.factor(year) + depth,
    time_varying = time_varying,
    spde = spde,
    time = time,
    family = tweedie(link = "log"),
    data = sub,
    anisotropy = TRUE,
    spatial_only = df$spatial_only[i],
    quadratic_roots = TRUE
  ), 
  silent=TRUE)
  
  if(class(m)!="try-error") saveRDS(m, file=paste0("output/wc/jacox_ROMS/model_ROMS_1_",i,".rds"))
  
  print(paste("species",i,df$covariate[i],sep=' ')) # just tracking the date of each haul to match to each var
}
Sys.time() - start.time

# Time difference of 2.901493 hours
