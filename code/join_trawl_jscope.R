# join empirical observations from trawl survey with environmental predictions from J-SCOPE

library(dplyr)

# read in pre-processed data sets
trawl_dat <- readRDS("survey_data/joined_nwfsc_data.rds")
jscope_dat <- readRDS("siedlecki_ROMS_output/joined_jscope_trawl_2009_2015.rds")

#define years of interest
years <- c(2009:2015)

# harmonize names, years, and drop corrupted trawl_id
trawl_dat$date <- gsub("-", "", as.POSIXlt(trawl_dat$date, format = "%Y-%b-%d"))

jscope_dat <- as.data.frame(jscope_dat) %>% 
  mutate(date = as.character(Date), year = as.integer(Year)) %>%
  select(year = Year, date, latitude_dd = Lat, longitude_dd = Lon,
         depth_model = "Model Depth", o2_model = "Model O2 (ml/l)", 
         sal_model = "Model Sal", temp_model = "Model Temp (°C)") %>%
  filter(year %in% years)

# sanity check on length of outputs from jscope relative to trawl survey observations
identical(nrow(jscope_dat), nrow(filter(trawl_dat, year %in% years, latitude_dd >= 43)))

# join
intersect_dat <- left_join(filter(trawl_dat, year %in% years, latitude_dd >= 43), jscope_dat) %>%
  filter(!is.na(o2_model))

# check for right number of distinct hauls (trawl_id) for a given species (as of 4/28/21, still incorrect!)
sp_dat <- filter(intersect_dat, scientific_name == "sebastes elongatus")
sum(unique(sp_dat$trawl_id))

saveRDS(intersect_dat, "survey_data/fully_joined_jscope_trawl_2009_2015.rds")
