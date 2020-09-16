library(ggplot2)
library(dplyr)

# read in estimates from all the ROMS models
region="wc" #c("goa","wc")[2]
df = readRDS(paste0("output/",region,"/jacox_ROMS/models_ROMS.RDS"))

for(i in 1:nrow(df)) {
  fname = paste0("output/",region,"/jacox_ROMS/model_ROMS_1_",i,".rds")
  if(file.exists(fname)) {
    
  m = readRDS(fname)
  sd_report <- summary(m$sd_report)
  params <- as.data.frame(sd_report[grep("quadratic", row.names(sd_report)), ])
  df$low[i] = params["quadratic_low","Estimate"]
  df$low_se[i] = params["quadratic_low","Std. Error"]
  df$hi[i] = params["quadratic_hi","Estimate"]
  df$hi_se[i] = params["quadratic_hi","Std. Error"]
  df$range[i] = params["quadratic_range","Estimate"]
  df$range_se[i] = params["quadratic_range","Std. Error"]
  df$peak[i] = params["quadratic_peak","Estimate"]
  df$peak_se[i] = params["quadratic_peak","Std. Error"]
  df$reduction[i] = params["quadratic_reduction","Estimate"]
  df$reduction_se[i] = params["quadratic_reduction","Std. Error"]  
  }
}

# remove junk columns and add a column descriminating ROMS vs in-situ
df_roms = dplyr::select(df, -spatial_only, -time_varying) %>% 
  dplyr::mutate("covariate"="ROMS")

# copy the same block but read in estimates from empirical
region="wc" #c("goa","wc")[2]
df = readRDS(paste0("output/",region,"/2003-2010 empirical/models_2003-2010.RDS"))
df$low = df$low_se = df$hi = df$hi_se = df$range = df$range_se = df$peak = df$peak_se = df$reduction = df$reduction_se
for(i in 1:nrow(df)) {
  fname = paste0("output/",region,"/2003-2010 empirical/model_",i,"_2003-2010.rds")
  if(file.exists(fname)) {
    
    m = readRDS(fname)
    sd_report <- summary(m$sd_report)
    params <- as.data.frame(sd_report[grep("quadratic", row.names(sd_report)), ])
    df$low[i] = params["quadratic_low","Estimate"]
    df$low_se[i] = params["quadratic_low","Std. Error"]
    df$hi[i] = params["quadratic_hi","Estimate"]
    df$hi_se[i] = params["quadratic_hi","Std. Error"]
    df$range[i] = params["quadratic_range","Estimate"]
    df$range_se[i] = params["quadratic_range","Std. Error"]
    df$peak[i] = params["quadratic_peak","Estimate"]
    df$peak_se[i] = params["quadratic_peak","Std. Error"]
    df$reduction[i] = params["quadratic_reduction","Estimate"]
    df$reduction_se[i] = params["quadratic_reduction","Std. Error"]  
  }
}

# remove junk columns and add a 
df = dplyr::select(df, -spatial_only, -time_varying) %>% 
  dplyr::mutate("covariate"="trawl survey")


# rbind the 2 datasets together
df = rbind(df, df_roms)

# filter out cases where models converged
df = df[complete.cases(df),]

# probably need to make decision about adding depth (or not). here
# we'll leave it in initially

pdf(paste0("plots/",region,"-temp_range2.pdf"))
level_order = dplyr::filter(df, !is.na(range), covariate=="trawl survey",
  depth_effect == TRUE, range_se < 1) %>%
  dplyr::arrange(range) %>% select(species)
dplyr::filter(df, !is.na(range), range_se < 1,depth_effect == TRUE,
  species %in% level_order$species) %>% 
  ggplot(aes(factor(species, level=level_order$species), range,col=covariate)) +
  geom_pointrange(aes(ymin=range-2*range_se, 
    ymax=range+2*range_se)) +
  coord_flip() + xlab("Species") + ylab("Range") + 
  ggtitle(paste0("Temperature range - ",region," survey"))
dev.off()


# also plot the empirical range of temps observed in both, 2003-2010
dat = readRDS("survey_data/joined_nwfsc_data.rds")
temp_survey = dat %>% dplyr::filter(species %in% level_order$species,
  year %in% seq(2003,2010)) %>% 
  group_by(species) %>% 
  summarize(temp_survey=diff(range(temp[which(cpue_kg_km2>0)],na.rm=T)))

ROMS.names <- c("ROMS_oxygen_bottom_era5_monthly","ROMS_temp_bottom_era5_monthly")
ROMS.RDS.names <- paste0("survey_data/joined_nwfsc_data",ROMS.names,".rds")

#dat = readRDS("survey_data/joined_nwfsc_data.rds")
dat_temp_bottom <- readRDS(ROMS.RDS.names[2])

temp_roms = dat_temp_bottom %>% dplyr::filter(species %in% level_order$species,
  year %in% seq(2003,2010)) %>% 
  group_by(species) %>% 
  summarize(temp_roms=diff(range(ROMS_temp_bottom_era5_monthly[which(cpue_kg_km2>0)],na.rm=T)))

temp_survey = dplyr::left_join(temp_survey, temp_roms)

pdf(paste0("plots/",region,"-survey_range_vs_roms_range.pdf"))
ggplot(temp_survey, aes(temp_survey, temp_roms,label=species)) + 
  geom_text() + geom_abline(intercept=0,slope=1,col="red")
dev.off()