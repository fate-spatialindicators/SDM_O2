#devtools::install_github("pbs-assess/sdmTMB")

source("code/mi_functions.R")
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)
library(rgdal)
library(future)

bc_dat <- load_bc_data(species.name = "sablefish", survey.name = "SYN WCVI")
wc_dat <- load_data()
bc_dat$Survey <- rep("BC", nrow(bc_dat))
wc_dat$Survey<- rep("WC", nrow(wc_dat))
dat <- full_join(bc_dat, wc_dat)
  
  # rescale variables
  dat$log_depth_scaled = scale(log(dat$depth))
  dat$log_depth_scaled2 = dat$log_depth_scaled^2
  mean.temp <- mean(dat$temp)
  sd.temp <- sd(dat$temp)
  mean.po2 <- mean(dat$po2)
  std.po2 <- sd(dat$po2)
  mean.do <- mean(dat$o2)
  std.do <- sd(dat$o2)
  mean.mi <- mean(dat$mi)
  std.mi <- sd(dat$mi)
  
  dat$temp = scale(dat$temp)
  dat$mi = scale(dat$mi)
  dat$o2 = scale(dat$o2)
  dat$po2 <- scale(dat$po2)
  dat$X <- dat$longitude
  dat$Y <- dat$latitude
  
  
  # run models for each combination of settings/covariates in df ------------
  
  use_cv = FALSE # specify whether to do cross validation or not
  use_AIC <- TRUE
  spde <- make_mesh(data = dat, xy_cols = c("X", "Y"), n_knots = 250) # choose # knots
  
  m_df <- get_models();
  AICmat <- dAIC <- tweedie_dens <- matrix(NA, nrow = length(m_df), ncol = 1) # set up array for tweedie density
  rownames(AICmat) <- rownames(dAIC) <- rownames(tweedie_dens) <- m_df
  plan(multisession)
  
  # fit models with typical approach
  for(i in 1:length(m_df)){
    print(paste0("model # ", i, " of ", length(m_df)))
    
    formula = paste0("cpue_kg_km2 ~ -1 + as.factor(Survey)+", m_df[i])
    
    # fit model with or without cross-validation
    if(use_cv==TRUE) {
      m <- sdmTMB_cv(
        formula = as.formula(formula),
        data = dat,
        time = NULL,
        spde = spde,
        k_folds = 5,
        family = tweedie(link = "log"),
        anisotropy = TRUE,
        spatial_only = TRUE,
        silent = TRUE)
      saveRDS(m, file = paste0("output/combined/cv/model_",i,"_MI_cv.rds"))
      tweedie_dens[i] = sum(m$fold_loglik)
      
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
        saveRDS(m, file = paste0("output/combined/model_",i,"_MI.rds"))
      }
      
    }
  }
  
  # If using AIC, calculate AIC and dAIC ------------------------------------
  if (use_AIC) {
    for (i in 1:length(m_df)) {
      filename <- paste0("output/combined/model_",i,"_MI.rds")
      m <- readRDS(filename)
      AICmat[i,1] <-AIC(m)
    }
  }
  
  
  
  
  dAIC[,1] <- AICmat[,1] - min(AICmat[,1])
  dAIC[,1] <-as.numeric(sprintf(dAIC,fmt = '%.2f'))
  dAIC
  