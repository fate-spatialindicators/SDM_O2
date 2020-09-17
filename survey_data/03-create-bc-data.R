library(dplyr)
library(ggplot2)

spp <- tolower(c(
  "Dover sole", "arrowtooth flounder", "rex sole",
  "English sole", "sablefish", "Pacific cod", "north pacific spiny dogfish",
  "longnose skate", "big skate", "Pacific ocean perch", "pacific halibut"
))

# Trawl survey data -------------------------------------------------------

.file <- "survey_data/bc-synoptic-tows.rds"
if (!file.exists(.file)) {
  dt <- list()
  for (i in seq_along(spp)) {
    .spp <- gsub(" ", "-", gsub("\\/", "-", tolower(spp[i])))
    # dt[[i]] <- gfdata::get_survey_sets(.spp, ssid = c(1, 3, 4, 16),
    #   join_sample_ids = TRUE)
    # or {
    dt[[i]] <- readRDS(paste0("../../gfs/report/data-cache/", .spp, ".rds"))
    dt[[i]] <- dt[[i]]$survey_sets
    # }
  }
  dt <- dplyr::bind_rows(dt)
  dt <- dplyr::filter(dt, survey_series_id %in% c(1, 3, 4, 16)) %>%
    dplyr::select(year,
      ssid = survey_series_id, 
      survey = survey_abbrev, species = species_common_name, longitude,
      latitude, depth_m, density_kgpm2, fishing_event_id
    )
  saveRDS(dt, file = .file)
} else {
  dt <- readRDS(.file)
}

# Environmental data ------------------------------------------------------

# New data:
# denv <- gfdata::get_sensor_data_trawl(ssid = c(1, 3, 4, 16), spread_attributes = FALSE)
# saveRDS(denv, "survey_data/bc-synoptic-env-raw.rds")
denv <- readRDS("survey_data/bc-synoptic-env-raw.rds")

# Old data:
# denv_old <- gfdata::get_table("FE_SALINITY_DO")
# saveRDS(denv_old, "survey_data/bc-do-salinity-table.rds")
denv_old <- readRDS("survey_data/bc-do-salinity-table.rds")
names(denv_old) <- tolower(names(denv_old))
denv_old <- denv_old %>%
  group_by(fishing_event_id) %>%
  summarise(
    avg = mean(do, na.rm = TRUE),
    min = min(do, na.rm = TRUE),
    max = max(do, na.rm = TRUE),
    start_time = min(fe_event_time),
    end_time = max(fe_event_time),
    year = lubridate::year(fe_event_time),
  .groups = "drop_last") %>%
  mutate(attribute = "do_mlpL") %>%
  distinct() 
denv_old

denv2 <- select(denv, fishing_event_id, attribute, avg, min, max, start_time, end_time, year)

denv2 <- bind_rows(denv2, denv_old)

denv2 %>% select(-start_time) %>% tidyr::pivot_wider(id_cols = c(year, fishing_event_id), names_from = attribute, values_from = avg) %>% head()

do <- denv2 %>% filter(attribute == "do_mlpL") %>%
  rename(do = avg) %>%
  select(year, fishing_event_id, do)

temp <- denv2 %>% filter(attribute == "temperature_C") %>%
  rename(temperature = avg) %>%
  select(year, fishing_event_id, temperature)

d <- left_join(dt, do)
d <- left_join(d, temp)

sum(is.na(d$do))
sum(is.na(d$temperature))

group_by(d, year) %>% summarise(frac_na = sum(is.na(do)) / n())
group_by(d, year) %>% summarise(frac_na = sum(is.na(temperature)) / n())

ggplot(d, aes(longitude, latitude, colour = do)) + geom_point() +
  facet_wrap(~year)

# ignore 2007, bad data:
d <- mutate(d, do = if_else(year != 2007, do, NA_real_))

ggplot(d, aes(longitude, latitude, colour = temperature)) + geom_point() +
  facet_wrap(~year)

env_data <- select(d, fishing_event_id, do, temperature, depth_m, longitude, latitude) %>% 
  mutate(fishing_event_id = as.integer(fishing_event_id)) %>% 
  distinct()
saveRDS(env_data, "survey_data/bc-synoptic-env.rds")

trawl_data <-  select(d, year, fishing_event_id, survey, species, density_kgpm2) %>%
  mutate(fishing_event_id = as.integer(fishing_event_id)) %>% 
  mutate(year = as.integer(year)) %>% 
  distinct()
saveRDS(trawl_data, "survey_data/bc-synoptic-trawls.rds")

# read and join with
env_data <- readRDS("survey_data/bc-synoptic-env.rds")
trawl_data <- readRDS("survey_data/bc-synoptic-trawls.rds")
dat <- left_join(trawl_data, env_data, by = "fishing_event_id")
