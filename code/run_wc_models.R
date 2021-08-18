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
libary(ggplot2)

# Set Run Specifications --------------------------------------------------
spc <- "sablefish"
fit.model <- F # do you want to fit SDM to environmental variables and impute missing values?
years <- 2010:2015 # designate years to use, must be 2010:2015 if using trawl-based data
constrain_latitude <- F # do you want to constraint trawl data N of 43 degrees latitude
compare_sources <- F # do you want to run models comparing trawl and J-SCOPE covariates?
no_depth <- FALSE # Do you want to run models w/out a depth effect?
use_cv = TRUE # specify whether to do cross validation or not
use_AIC = FALSE # specify whether to use AIC
use_jscope <- F # specify whether to only use J-SCOPE based estimates.  Overrides compare_sources and fit.model

# parallel cross-validation:
if (use_cv) {
  is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
  is_unix <- .Platform$OS.type == "unix"
  if (!is_rstudio && is_unix) plan(multicore, workers = 6L) else plan(multisession)
}

# load and scale data -----------------------------------------------------
if (!use_jscope) {
  dat <- load_data(fit.model, spc, constrain_latitude)
  if (compare_sources) {
    # Load J-SCOPE based, and merge into single data set
    jscope_dat <- load_data_jscope(spc, years)
    dat <- left_join(x = dat, y = jscope_dat, by = "trawl_id")
    # remove missing tows
    dat <-
      dplyr::filter(dat,!is.na(cpue_kg_km2.y),!is.na(cpue_kg_km2.x))
    # simplify and rename variables
    dat <- dplyr::select(
      dat,
      year = year.x,
      cpue_kg_km2 = cpue_kg_km2.x,
      julian_day = julian_day.x,
      longitude = longitude.x,
      latitude = latitude.x,
      o2_trawl = o2.x,
      temp_trawl = temp.x,
      po2_trawl = po2.x,
      mi_trawl = mi.x,
      o2_jscope = o2.y,
      temp_jscope = temp.y,
      po2_jscope = po2.y,
      mi_jscope = mi.y,
      depth = depth.x
    )
    # rescale variables
    dat$temp_trawl <- as.numeric(scale(dat$temp_trawl))
    dat$mi_trawl <- as.numeric(scale(dat$mi_trawl))
    dat$po2_trawl <- as.numeric(scale(dat$po2_trawl))
    dat$temp_jscope <- as.numeric(scale(dat$temp_jscope))
    dat$mi_jscope <- as.numeric(scale(dat$mi_jscope))
    dat$po2_jscope <- as.numeric(scale(dat$po2_jscope))
    dat$log_depth_scaled <- scale(log(dat$depth))
    dat$log_depth_scaled2 <- with(dat, log_depth_scaled ^ 2)
    dat$log_depth_scaled3 <- with(dat, log_depth_scaled^3)
    dat$jday_scaled <- scale(dat$julian_day)
    dat$jday_scaled2 <- with(dat, jday_scaled ^ 2)
    dat$X <- dat$longitude
    dat$Y <- dat$latitude
  }
  if (!compare_sources) {
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
  }
}

if (use_jscope) {
  dat <- load_data_jscope(spc, years)
  dat$temp <- as.numeric(scale(dat$temp))
  dat$mi <- as.numeric(scale(dat$mi))
  dat$po2 <- as.numeric(scale(dat$po2))
  dat$log_depth_scaled <- scale(log(dat$depth))
  dat$log_depth_scaled2 <- with(dat, log_depth_scaled ^ 2)
  dat$log_depth_scaled3 <- with(dat, log_depth_scaled^3)
  dat$jday_scaled <- scale(dat$julian_day)
  dat$jday_scaled2 <- with(dat, jday_scaled ^ 2)
  dat$X <- dat$longitude
  dat$Y <- dat$latitude
}
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

if (use_cv) {
  tweedie_dens <- matrix(NA_real_, nrow = length(m_df), ncol = 3L)
} else {
  tweedie_dens <- matrix(NA_real_, nrow = length(m_df), ncol = 1L)
}

# fit models and save files to output/wc folder
for (j in seq_len(ncol(tweedie_dens))) {
  set.seed(j * 102849)
  dat$fold_ids <- sample(seq_len(12), nrow(dat), replace = TRUE)
  for (i in seq(1, length(m_df))) {
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
    if (use_cv) {
      m <- sdmTMB_cv(
        formula = as.formula(formula),
        data = dat,
        time = NULL,
        spde = spde,
        fold_ids = dat$fold_ids,
        family = tweedie(link = "log"),
        anisotropy = TRUE,
        spatial_only = TRUE,
        silent = TRUE,
        control = sdmTMBcontrol(nlminb_loops = 1, newton_loops = 1)
      )
      dir.create("output/wc/cv/", showWarnings = FALSE)
      saveRDS(m, file = paste0("output/wc/cv/model_", i, j, "_", spc, "_cv.rds"))
      tweedie_dens[i,j] <- sum(m$fold_loglik)
      
    } else {
      m <- try(sdmTMB(
        formula = as.formula(formula),
        data = dat,
        time = NULL,
        reml = FALSE,
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
          saveRDS(m, file = paste0("output/wc/model_", i, "_", spc, ".rds"))
        if(!no_depth&use_jscope) saveRDS(m, file = paste0("output/wc/model_", i, "_", spc, "_jscope.rds"))
      }
    }
  }
}

# If using AIC, calculate AIC and dAIC ------------------------------------
if (use_AIC) {
  for (i in 1:length(m_df)) {
    if (no_depth)
      filename <- paste0("output/wc/model_", i, "_", spc, "_nodepth.rds")
    if (!no_depth&!use_jscope)
      filename <- paste0("output/wc/model_", i, "_", spc, ".rds")
    if(!no_depth&use_jscope)
      filename <- paste0("output/wc/model_", i, "_", spc, "_jscope.rds")
    m <- readRDS(filename)
    AICmat[i, 1] <- AIC(m)
  }
  
  # calculate and print out delta AIC table
  dAIC[, 1] <- AICmat[, 1] - min(AICmat[, 1])
  dAIC[, 1] <- as.numeric(sprintf(dAIC, fmt = '%.2f'))
  print(spc)
  dAIC
  
  # print best model
  best.index <- which(dAIC==0)
}

if (use_cv) saveRDS(tweedie_dens, file = "output/wc/cv/tweedie_dens.rds")

if (!use_cv) {
  if (no_depth)
    filename <- paste0("output/wc/model_", best.index, "_", spc, "_nodepth.rds")
  if (!no_depth&!use_jscope)
    filename <- paste0("output/wc/model_", best.index, "_", spc, ".rds")
  if(!no_depth&use_jscope) 
    filename <- paste0("output/wc/model_", best.index, "_", spc, "_jscope.rds")
  
  m <- readRDS(filename)
  summary(m)
}

if (use_cv) {
  tweedie_dens <- readRDS("output/wc/cv/tweedie_dens.rds")
  xx <- data.frame(f = m_df, Var1 = seq_along(m_df), stringsAsFactors = FALSE)
  reshape2::melt(tweedie_dens) %>%
    left_join(xx) %>% 
    rename(split = Var2) %>% 
    group_by(split) %>% 
    mutate(value = value - max(value)) %>%
    ggplot(aes(f, value, colour = as.factor(split), group = split)) + 
    geom_point() +
    geom_line() +
    coord_flip() 
}