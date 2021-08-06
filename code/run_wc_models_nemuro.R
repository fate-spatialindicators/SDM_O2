# Master file to fit SDM to west coast bottom trawl survey, includes options to
# (1) use Trawl-based or J-SCOPE based estimates
# (2) Impute missing trawl-based oxygen and salinity measurements
# (3) Do model selection comparing trawl-based or J-SCOPE based environmental covariates


# Install latest version of sdmTMB ----------------------------------------
#devtools::install_github("pbs-assess/sdmTMB")
# Initialize with packages and functions ------------------------------------------------
rm(list = ls())
source("code/mi_functions.R")
library(sdmTMB)
library(raster)
library(dplyr)
library(sp)
library(gsw)
library(rgdal)
library(future)



# Set Run Specifications --------------------------------------------------
spc <- "sablefish"
fit.model <- F # do you want to fit SDM to environmental variables and impute missing values?
constrain_latitude <- F # do you want to constraint trawl data N of 43 degrees latitude
compare_sources <- F # do you want to run models comparing trawl and J-SCOPE covariates?
no_depth <- FALSE # Do you want to run models w/out a depth effect?
use_cv = FALSE # specify whether to do cross validation or not
use_AIC = TRUE # specify whether to use AIC
use_jscope <- F # specify whether to only use J-SCOPE based estimates.  Overrides compare_sources and fit.model


# load and scale data -----------------------------------------------------

  dat <- load_data_nemuro(fit.model, spc, constrain_latitude)
  
  
    dat$temp <- (scale(dat$temp))
    dat$mi <- (scale(dat$mi))
    dat$po2 <- (scale(dat$po2))
    dat$log_depth_scaled <- scale(log(dat$depth))
    dat$log_depth_scaled2 <- with(dat, log_depth_scaled ^ 2)
    dat$log_depth_scaled3 <- with(dat, log_depth_scaled^3)
    dat$jday_scaled <- scale(dat$julian_day)
    dat$jday_scaled2 <- with(dat, jday_scaled ^ 2)
    dat$X <- dat$longitude
    dat$Y <- dat$latitude


# make year a factor
dat$year <- as.factor(dat$year)

# Run alternative models  -----------------------------------------------------
spde <-
  make_mesh(data = dat,
            xy_cols = c("X", "Y"),
            n_knots = 250) # choose # knots

# get list of models.  if use_scope is T, then always use "get_models()".
if (compare_sources & !use_jscope)
  m_df <- get_models_compare()
if (!compare_sources | use_jscope)
  m_df <- get_models()

AICmat <-
  dAIC <-
  matrix(NA, nrow = length(m_df), ncol = 1) # set up array for output
rownames(AICmat) <- rownames(dAIC) <- m_df


# fit models and save files to output/wc folder
for (i in 1:length(m_df)) {
  print(paste0("model # ", i, " of ", length(m_df)))
  formula = paste0("cpue_kg_km2 ~ -1 +", m_df[i])
  if (no_depth)
    formula <-
    gsub(
      pattern = "log_depth_scaled + log_depth_scaled2 +",
      x = formula,
      replacement = "",
      fixed = T
    )
  # fit model with or without cross-validation
  if (use_cv == TRUE) {
    m <- sdmTMB_cv(
      formula = as.formula(formula),
      data = dat,
      time = NULL,
      spde = spde,
      k_folds = 5,
      family = tweedie(link = "log"),
      anisotropy = TRUE,
      spatial_only = TRUE,
      silent = TRUE
    )
    saveRDS(m, file = paste0("output/wc/cv/model_", i, "_", spc, "_cv.rds"))
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
    ),
    silent = TRUE)
    
    if (class(m) != "try-error") {
      if (no_depth)
        saveRDS(m, file = paste0("output/wc/model_", i, "_", spc, "_nodepth.rds"))
      if (!no_depth&!use_jscope)
        saveRDS(m, file = paste0("output/wc/model_", i, "_", spc, "_nemuro.rds"))
      if(!no_depth&use_jscope) saveRDS(m, file = paste0("output/wc/model_", i, "_", spc, "_jscope.rds"))
    }
  }
}

# If using AIC, calculate AIC and dAIC ------------------------------------
if (use_AIC) {
  for (i in 1:length(m_df)) {
    if (no_depth)
      filename <- paste0("output/wc/model_", i, "_", spc, "_nodepth.rds")
    if (!no_depth&!use_jscope)
      filename <- paste0("output/wc/model_", i, "_", spc, "_nemuro.rds")
    if(!no_depth&use_jscope)
      filename <- paste0("output/wc/model_", i, "_", spc, "_jscope.rds")
    m <- readRDS(filename)
    AICmat[i, 1] <- AIC(m)
  }
}

# calculate and print out delta AIC table

dAIC[, 1] <- AICmat[, 1] - min(AICmat[, 1])
dAIC[, 1] <- as.numeric(sprintf(dAIC, fmt = '%.2f'))
print(spc)
dAIC

# print best model
best.index <- which(dAIC==0)
if (no_depth)
  filename <- paste0("output/wc/model_", best.index, "_", spc, "_nodepth.rds")
if (!no_depth&!use_jscope)
  filename <- paste0("output/wc/model_", best.index, "_", spc, "_nemuro.rds")
if(!no_depth&use_jscope) 
  filename <- paste0("output/wc/model_", best.index, "_", spc, "_jscope.rds")
  
m <- readRDS(filename)
summary(m)

