# join empirical observations from trawl survey with environmental predictions from J-SCOPE

library(dplyr)
library(nwfscSurvey)

# read in pre-processed data sets
haul <- nwfscSurvey::PullHaul.fn(SurveyName = "NWFSC.Combo")
haul <- plyr::rename(haul, replace=c("salinity_at_gear_psu_der" = "sal", 
                                     "temperature_at_gear_c_der" = "temp", 
                                     "o2_at_gear_ml_per_l_der" = "o2",
                                     "depth_hi_prec_m" = "depth"))
haul$year = as.numeric(substr(haul$date_yyyymmdd,1,4))
haul$month = as.numeric(substr(haul$date_yyyymmdd,5,6))
haul$day = as.numeric(substr(haul$date_yyyymmdd,7,8))
#define years of interest
years <- c(2003:2018)

trawl_dat <- haul %>%
  filter(year %in% years, latitude_dd >= 43)

# make tow_id 
tow_id <- with(trawl_dat, rank(trawl_id))
trawl_dat$tow_id <- tow_id

#trawl_dat <- readRDS("survey_data/joined_nwfsc_data.rds")
jscope_dat <- readRDS("siedlecki_ROMS_output/jscope_estimates.rds")
jscope_dat$year <- as.numeric(substr(jscope_dat$date,1,4))
jscope_dat$month <- as.numeric(substr(jscope_dat$date,5,6))
jscope_dat$day <- as.numeric(substr(jscope_dat$date,7,8))


# harmonize names, years, and drop corrupted trawl_id
jscope_dat <- as.data.frame(jscope_dat) %>% 
  dplyr::select(year, month, day, 
         depth_model, o2_model, 
         sal_model, temp_model, pass, tow_id) %>%
  filter(year %in% years)


intersect_dat <- (inner_join(trawl_dat, jscope_dat, by = c("tow_id")))

intersect_dat <- dplyr::select(intersect_dat, year = year.x, month = month.x, day = day.x, date_yyyymmdd, longitude_dd, latitude_dd,
                               o2, temp, sal, o2_model, temp_model, sal_model, depth, pass=pass.x, performance, trawl_id)


saveRDS(intersect_dat, "survey_data/joined_jscope_trawl_2003_2018.rds")


