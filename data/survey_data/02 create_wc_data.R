library(nwfscSurvey)
library(sdmTMB)
library(dplyr)
library(stringr)

# bring in common names
#UrlText <- "https://www.nwfsc.noaa.gov/data/api/v1/source/trawl.catch_fact/selection.json?filters=project=Groundfish%20Slope%20and%20Shelf%20Combination%20Survey,date_dim$year>=2003&variables=common_name,scientific_name"
#DataPull <- try(jsonlite::fromJSON(UrlText))
#spec_names = group_by(DataPull, common_name) %>% 
#  dplyr::summarize(scientific_name = scientific_name[1])
#spec_names$scientific_name = tolower(spec_names$scientific_name)
#saveRDS(spec_names,"nwfsc_lookup.rds")

spec_names = readRDS("survey_data/nwfsc_species_lookup.rds")

catch = readRDS("survey_data/wcbts_catch_2019-08-01.rds")
names(catch) = tolower(names(catch))
catch$date = as.character(catch$date)
catch$trawl_id = as.numeric(catch$trawl_id)
catch$scientific_name = tolower(catch$scientific_name)
catch = dplyr::left_join(catch, spec_names)
catch$common_name = tolower(catch$common_name)

haul = readRDS("survey_data/wcbts_haul_2019-08-01.rds")

#haul = dplyr::select(haul, date_yyyymmdd,depth_hi_prec_m,latitude_dd,longitude_dd,temperature_at_gear_c_der)
#write.table(haul, "data_for_mike.csv",sep=",",row.names=FALSE)
haul$year = as.numeric(substr(haul$date_yyyymmdd,1,4))
haul$month = as.numeric(substr(haul$date_yyyymmdd,5,6))
haul$day = as.numeric(substr(haul$date_yyyymmdd,7,8))

haul = dplyr::rename(haul, 
  o2 = o2_at_gear_ml_per_l_der,
  degc = temperature_at_gear_c_der,
  sal = salinity_at_gear_psu_der,
  depthm=depth_hi_prec_m) %>% 
  dplyr::select(o2,degc,sal,depthm,latitude_dd,longitude_dd,
    performance,trawl_id)

dat = dplyr::left_join(catch, haul)

# WC names in cope and haltuch 2012
cope_haltuch = c("aurora rockfish", "big skate", "bigfin eelpout",
  "black eelpout", "brown cat shark", "california slickhead",
  "canary rockfish", "chilipepper", "darkblotched rockfish",
  "deepsea sole", "dover sole", "english sole", "giant grenadier",
  "greenstriped rockfish", "halfbanded rockfish", "lingcod", 
  "longnose skate", "longspine thornyhead", "butterfish unident.",
  "pacific flatnose", "pacific grenadier", "pacific hake",
  "pacific sanddab", "petrale sole", "pink seaperch",
  "rex sole", "sablefish", "sandpaper skate", "sharpchin rockfish",
  "shortbelly rockfish", "shortspine thornyhead", "slender sole",
  "pacific spiny dogfish", "splitnose rockfish", "spotted ratfish",
  "stripetail rockfish", "white croaker", "yellowtail rockfish")

# these are additional species on the prioritization spreadsheet
fram = c(cope_haltuch, "vermilion and sunset rockfish",
  "black rockfish", "cowcod","copper rockfish",
  "brown rockfish","quillback rockfish",
  "redbanded rockfish","tree rockfish",
  "squarespot rockfish","starry rockfish",
  "speckled rockfish","rougheye and blackspotted rockfish",
  "shortraker rockfish","flathead sole",
  "widow rockfish","kelp greenling",
  "olive rockfish","blue rockfish",
  "kelp rockfish","cabezon",
  "sand sole","flag rockfish",
  "starry flounder","rock sole unident.",
  "greenspotted rockfish",
  "honeycomb rockfish",
  "California scorpionfish",
  "blackgill rockfish")
# filter by 'well sampled wc species'
#dat = dplyr::filter(dat, common_name %in% cope_haltuch)

# filter by species that occur in 10% of hauls
threshold = 0.1

keep = dat %>% 
  mutate(occur = ifelse(cpue_kg_km2 > 0,1,0)) %>%
  group_by(year, common_name) %>% 
  summarize(p = sum(occur)/n()) %>% 
  group_by(common_name) %>% 
  summarize(mean_p = mean(p, na.rm=T)) %>% 
  filter(mean_p >= threshold)

keep = c(keep$common_name, 
  fram[which(fram %in% keep$common_name==FALSE)])

# remove urchins, stars, etc
spp_to_keep = c("urchin|star|anemone|cucumber|sea pen|salps|sponge|snail|shab|slug|jellyfish|squid|shrimp|hagfish|pasiphaeid|smelt|tongue|tritonia")
spp_to_keep = keep[-grep(spp_to_keep,keep)]
dat = dplyr::filter(dat, common_name %in% spp_to_keep)

dat = dplyr::rename(dat, 
  temp=degc,depth=depthm,species=common_name) %>% 
  select(-depth_m)
dat$survey = "nwfsc.combo"

saveRDS(dat, "survey_data/joined_nwfsc_data.rds")
