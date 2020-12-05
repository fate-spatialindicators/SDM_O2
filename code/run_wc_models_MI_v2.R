# Install latest version of sdmTMB ----------------------------------------
devtools::install_github("pbs-assess/sdmTMB")

# Initialize with packages and functions ------------------------------------------------
rm(list = ls())
source("code/mi_functions.R")
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)
library(rgdal)

# load and scale data -----------------------------------------------------

dat <- load_data()
dat$log_depth_scaled <- scale(log(dat$depth))
dat$log_depth_scaled2 <- dat$log_depth_scaled^2
dat$jday_scaled <- scale(dat$julian_day)
dat$jday_scaled2 <- dat$jday_scaled^2
dat$X <- dat$longitude
dat$Y <- dat$latitude

# rescale variables
dat$mi <- as.numeric(scale(dat$mi))
dat$o2 <- as.numeric(scale(dat$o2))
dat$po2 <- as.numeric(scale(dat$po2))


# Run alternative models  -----------------------------------------------------

use_cv = FALSE # specify whether to do cross validation or not
use_AIC = TRUE # specify whether to use AIC
spde <- make_mesh(data = dat, xy_cols = c("X", "Y"), n_knots = 250) # choose # knots

m_df <- get_models();
AICmat <- dAIC <- tweedie_dens <- matrix(NA, nrow = length(m_df), ncol = 1) # set up array for tweedie density
rownames(AICmat) <- rownames(dAIC) <- rownames(tweedie_dens) <- m_df


# fit models with typical approach
for(i in 1:length(m_df)){
  print(paste0("model # ", i, " of ", length(m_df)))
  
  formula = paste0("cpue_kg_km2 ~ 0 +", m_df[i])
  
  # fit model with or without cross-validation
  if(use_cv==TRUE) {
    m <- try(sdmTMB_cv(
      formula = as.formula(formula),
      data = dat,
      time = NULL,
      spde = spde,
      k_folds = 10,
      family = tweedie(link = "log"),
      anisotropy = TRUE,
      spatial_only = TRUE
    ), silent = TRUE)
    
    if(class(m)!="try-error") {
      saveRDS(m, file = paste0("output/wc/model_",i,"_MI_cv.rds"))
      tweedie_dens[i] = m$sum_loglik
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

# If using AIC, calculate AIC and dAIC ------------------------------------
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
