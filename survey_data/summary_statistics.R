# code to summarize mean lat, lon, depth, etc by species
library(dplyr)

dat = readRDS("survey_data/joined_nwfsc_data.rds")

g = dplyr::group_by(dat, species) %>% 
  summarize(mean_depth = mean(depth[which(cpue_kg_km2>0)],na.rm=T),
    min_depth = min(depth[which(cpue_kg_km2>0)],na.rm=T),
    max_depth = max(depth[which(cpue_kg_km2>0)],na.rm=T),
    lower95_depth = quantile(depth[which(cpue_kg_km2>0)],0.025),
    upper95_depth = quantile(depth[which(cpue_kg_km2>0)],0.975), 
    range_depth = max_depth-min_depth,
    range95_depth = upper95_depth - lower95_depth,
    mean_lat = mean(latitude_dd[which(cpue_kg_km2>0)],na.rm=T),
    mean_lon = mean(longitude_dd[which(cpue_kg_km2>0)],na.rm=T),
    weighted_depth = sum(depth*cpue_kg_km2,na.rm=T)/sum(cpue_kg_km2,na.rm=T),
    weighted_lat = sum(latitude_dd*cpue_kg_km2,na.rm=T)/sum(cpue_kg_km2,na.rm=T),
    weighted_lon = sum(longitude_dd*cpue_kg_km2,na.rm=T)/sum(cpue_kg_km2,na.rm=T))

d = dplyr::group_by(dat, species, year) %>% 
  summarize(p = length(which(cpue_kg_km2>0))/n(),
    n = length(which(cpue_kg_km2>0))) %>%
  ungroup() %>% 
  group_by(species) %>% 
  summarize(n_year = length(which(p > 0)), p = mean(p),
    min_n = min(n), mean_n = mean(n))

g = dplyr::left_join(g, d)
write.csv(as.data.frame(g), file="output/summary_statistics_wc.csv")

dat = readRDS("survey_data/joined_goa_data.rds")

g = dplyr::group_by(dat, species) %>% 
  summarize(mean_depth = mean(depth[which(cpue_kg_km2>0)],na.rm=T),
    min_depth = min(depth[which(cpue_kg_km2>0)],na.rm=T),
    max_depth = max(depth[which(cpue_kg_km2>0)],na.rm=T),
    lower95_depth = quantile(depth[which(cpue_kg_km2>0)],0.025),
    upper95_depth = quantile(depth[which(cpue_kg_km2>0)],0.975), 
    range_depth = max_depth-min_depth,
    range95_depth = upper95_depth - lower95_depth,
    mean_lat = mean(latitude_dd[which(cpue_kg_km2>0)],na.rm=T),
    mean_lon = mean(longitude_dd[which(cpue_kg_km2>0)],na.rm=T),
    weighted_depth = sum(depth*cpue_kg_km2,na.rm=T)/sum(cpue_kg_km2,na.rm=T),
    weighted_lat = sum(latitude_dd*cpue_kg_km2,na.rm=T)/sum(cpue_kg_km2,na.rm=T),
    weighted_lon = sum(longitude_dd*cpue_kg_km2,na.rm=T)/sum(cpue_kg_km2,na.rm=T))

d = dplyr::group_by(dat, species, year) %>% 
  summarize(p = length(which(cpue_kg_km2>0))/n(),
    n = length(which(cpue_kg_km2>0))) %>%
  ungroup() %>% 
  group_by(species) %>% 
  summarize(n_year = length(which(p > 0)), p = mean(p),
    min_n = min(n), mean_n = mean(n))

g = dplyr::left_join(g, d)
write.csv(as.data.frame(g), file="output/summary_statistics_goa.csv")
