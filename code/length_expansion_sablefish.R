## Length Expansion ##

library(sp)
library(ggplot2)
library(broom)
library(dplyr)
library(tidyr)
library(purrr)
library(KernSmooth)

# Species of interest and max. juvenile lengths (define ontogenetic classes)
species = read.csv("survey_data/species_list_revised.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name,
                        juv_threshold = max.length.cm)

# load, clean, and join data
bio = readRDS("data/survey_data/wcbts_bio_2019-08-01.rds")
haul = readRDS("data/survey_data/wcbts_haul_2019-08-01.rds")
catch = readRDS("data/survey_data/wcbts_catch_2019-08-01.rds")
names(catch) = tolower(names(catch))
names(bio) = tolower(names(bio))
names(haul) = tolower(names(haul))

bio$trawl_id = as.character(bio$trawl_id)
haul$trawl_id = as.character(haul$trawl_id)
haul$date_yyyymmdd = as.numeric(haul$date_yyyymmdd)
haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)

dat = dplyr::left_join(catch[,c("trawl_id","scientific_name","year","subsample_count",
                                "subsample_wt_kg","total_catch_numbers","total_catch_wt_kg","cpue_kg_km2")], haul) %>%
  dplyr::left_join(filter(bio, !is.na(length_cm))) %>%
  filter(performance == "Satisfactory")  %>%
  mutate(depth_m = depth_hi_prec_m)


  sci_name = "Anoplopoma fimbria"
  # filter out species of interest from joined (catch/haul/bio) dataset
  dat_sub = dplyr::filter(dat, scientific_name == sci_name)
  
  # fit length-weight regression by year to predict fish weights that have lengths only.
  # note a rank-deficiency warning may indicate there is insufficient data for some year/sex combinations (likely for unsexed group)

    fitted = dat_sub %>%
    filter(!is.na(length_cm), !is.na(weight_kg)) %>%
    dplyr::select(trawl_id,year,
           subsample_wt_kg, total_catch_wt_kg, area_swept_ha_der, cpue_kg_km2,
           individual_tracking_id, sex, length_cm, weight_kg) %>%
    group_nest(year)  %>%
    mutate(
      model = map(data, ~ lm(log(weight_kg) ~ log(length_cm), data = .x)),
      tidied = map(model, tidy),
      augmented = map(model, augment),
      predictions = map2(data, model, modelr::add_predictions)
    )
  
  # replace missing weights with predicted weights
  dat_pos = fitted %>%
    unnest(predictions) %>%
    dplyr::select(-data, -model, -tidied, -augmented) %>%
    mutate(weight = ifelse(is.na(weight_kg), exp(pred), weight_kg))
  
  trawlids <- unique(dat_pos$trawl_id)
  p <- data.frame(trawl_id = trawlids,
                             p1 = 0,
                             p2 = 0,
                             p3 = 0,
                             p4 = 0,
                            p1n = 0,
                            p2n = 0,
                            p3n = 0,
                            p4n = 0)
  
  sizethresholds <- c(0.5, 2, 6, 11)
  for (i in 1:length(trawlids)) {
  haul_sample<- dplyr::filter(dat_pos, trawl_id == trawlids[i])
  if(nrow(haul_sample) > 0 | var(haul_sample$weight >0)) {
  # fit kernel density to weight frequency
  smoothed_w <- bkde(haul_sample$weight, range.x = c(min(dat_pos$weight), max(dat_pos$weight)), bandwidth = 2)
  # make sure smoother predicts positive or zero density
  smoothed_w$y[smoothed_w$y<0] <- 0
  # calculate proportion by biomass and by number
  p_w_byweight <- smoothed_w$y * smoothed_w$x / sum(smoothed_w$x*smoothed_w$y)
  p_w_bynum <- smoothed_w$y / sum(smoothed_w$y)
  
  #p_w_byweight[p_w_byweight<0] <- 0
  #p_w_bynum[p_w_bynum<0] <- 0
  
  p1 <- sum(p_w_byweight[smoothed_w$x<=0.5])
  p2 <- sum(p_w_byweight[smoothed_w$x>0.5 & smoothed_w$x <=2])
  p3 <- sum(p_w_byweight[smoothed_w$x>2 & smoothed_w$x <=6])
  p4 <- sum(p_w_byweight[smoothed_w$x>6])
  
  p1n <- sum(p_w_bynum[smoothed_w$x<=0.5])
  p2n <- sum(p_w_bynum[smoothed_w$x>0.5 & smoothed_w$x <=2])
  p3n <- sum(p_w_bynum[smoothed_w$x>2 & smoothed_w$x <=6])
  p4n <- sum(p_w_bynum[smoothed_w$x>6])
  
  
  p[i,2:9] <- c(p1, p2, p3, p4, p1n, p2n, p3n, p4n)
  
  }
  else {
    indx <- which(sizethresholds>haul_sample$weight)
    p[i, min(indx)+1] <- 1
    p[i, min(indx)+5] <- 1
  }
  }
    
  
  # add hauls with zero catch back in
  absent = filter(dat_sub, cpue_kg_km2 == 0)
  trawlids <- unique(absent$trawl_id)
  absent.df <- data.frame(trawl_id = trawlids,
                          p1 = 0,
                          p2 = 0,
                          p3 = 0,
                          p4 = 0,
                          p1n = 0,
                          p2n = 0,
                          p3n = 0,
                          p4n = 0)
  
  all_hauls <- rbind(p, absent.df)
  all_hauls$trawl_id <- as.numeric(all_hauls$trawl_id)
  saveRDS(all_hauls, "data/survey_data/sablefish_size_dist.rds")
  
  
  