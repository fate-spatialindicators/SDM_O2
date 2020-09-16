#devtools::install_github("pbs-assess/sdmTMB")
rm(list = ls())
source("code/mi_functions.R")

library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)
library(ggplot2)
dat <- load_data()

# rescale variables
dat$depth = scale(log(dat$depth))
dat$temp = scale(dat$temp)

mean.po2 <- mean(dat$po2)
std.po2 <- sd(dat$po2)
mean.do <- mean(dat$o2)
std.do <- sd(dat$o2)
mean.mi <- mean(dat$mi)
std.mi <- sd(dat$mi)
dat$mi = scale(dat$mi)
dat$o2 = scale(dat$o2)
dat$po2 <- scale(dat$po2)


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


# run models for each combination of settings/covariates in df ------------

use_cv = FALSE # specify whether to do cross validation or not
use_AIC = TRUE # specify whether to use AIC
spde <- make_spde(x = dat$longitude, y = dat$latitude, n_knots = 250) # choose # knots

m_df <- get_models();
AICmat <- dAIC <- tweedie_dens <- matrix(NA, nrow = length(m_df), ncol =1 ) # set up array for tweedie density
rownames(AICmat) <- rownames(dAIC) <- rownames(tweedie_dens) <- m_df


# fit models with typical approach
for(i in 1:length(m_df)){
  print(paste0("model # ", i, " of ", length(m_df)))
  
  formula = paste0("cpue_kg_km2 ~ 0 +", m_df[i])
  
  # fit model with or without cross-validation
  if(use_cv==TRUE) {
    m <- try(sdmTMB_cv(
      formula = as.formula(formula),
      data = sub,
      x = "longitude", 
      y = "latitude",
      time = NULL,
      k_folds = 4,
      n_knots = 250,
      seed = 10,
      family = tweedie(link = "log"),
      anisotropy = TRUE,
      spatial_only = TRUE
    ), silent = TRUE)
    
    if(class(m)!="try-error") {
      saveRDS(m, file = paste0("output/wc/model_",i,"_MI_cv.rds"))
      m_df$tweedie_dens[i] = m$sum_loglik
    }
    
  } else {
    m <- try(sdmTMB(
      formula = as.formula(formula),
      data = dat,
      time = NULL,
      reml = TRUE,
      spde = spde,
      family = tweedie(link = "log"),
      anisotropy = TRUE,
      spatial_only = TRUE
    ), silent = TRUE)
    
    if(class(m)!="try-error") {
      saveRDS(m, file = paste0("output/wc/model_",i,"_MI.rds"))
    }
    
  }
}

# Use this if using AIC
if (use_AIC) {
for (i in 1:length(m_df)) {
  filename <- paste0("output/wc/model_",i,"_MI.rds")
  m <- readRDS(filename)
  AICmat[i,1] <-AIC(m)
}

dAIC[,1] <- AICmat[,1] - min(AICmat[,1])
dAIC[,1] <-as.numeric(sprintf(dAIC,fmt = '%.2f'))
dAIC
}
