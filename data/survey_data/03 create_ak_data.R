library(dplyr)

d = readRDS("survey_data/AK_BTS_all_spp.rds")

names(d) = tolower(names(d))

# rename to match nwfsc
d = dplyr::rename(d, species = common_name,
  scientific_name = species_name,
  latitude_dd = latitude,
  longitude_dd = longitude,
  cpue_kg_km2=cpue,
  temp=gear_temperature,
  depth=bottom_depth)

# filter out goa
d = dplyr::filter(d,survey=="GOA")

# filter out 0 depths
d = dplyr::filter(d,depth>0)

# filter out years < 1990 (gear changes) and 2001 (limited sampling)
d = dplyr::filter(d, year%in%c(1984,1987,2001)==FALSE)

saveRDS(d, "survey_data/joined_goa_data.rds")

