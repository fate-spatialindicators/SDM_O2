# Load functions and packages
#devtools::install_github("pbs-assess/sdmTMB")
source("code/mi_functions.R")
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)
library(ggplot2)
library(viridis)
library(gridExtra)
library(raster)
library(rasterize)
library(rgdal)

run_models <- T
fit_env <- F

# load models and data
formulas <- get_models()
dat <- load_data()


dat <- load_bc_data(species.name = "sablefish", survey.name = "SYN WCVI")

years <- unique(dat$year)


load("survey_data/synoptic_grid.rda")
wcvi_grid <- dplyr::filter(synoptic_grid, survey == "SYN WCVI")


dat$log_depth_scaled <- scale(log(dat$depth))
dat$log_depth_scaled2 <- dat$log_depth_scaled ^ 2
dat$mi_unscaled <- as.numeric(dat$mi)
dat$po2_unscaled <- as.numeric(dat$po2)
dat$o2_unscaled <- as.numeric(dat$o2)
dat$temp_unscaled <- as.numeric(dat$temp)

dat$mi <- as.numeric(scale(dat$mi))
dat$o2 <- as.numeric(scale(dat$o2))
dat$po2 <- as.numeric(scale(dat$po2))
dat$temp <- as.numeric(scale(dat$temp))
dat$X <- dat$longitude
dat$Y <- dat$latitude


c_spde <-
  make_mesh(data = dat,
            xy_cols = c("X", "Y"),
            n_knots = 250) # choose # knots


if (run_models) {
  if (fit_env) {
    # Make grid of environmental conditions - fit scaled coefficients
    temp_model <-
      sdmTMB(
        formula = temp ~ -1 + log_depth_scaled + log_depth_scaled2
        +  as.factor(year),
        data = dat,
        time = "year",
        spde = c_spde,
        anisotropy = TRUE,
        silent = TRUE,
        spatial_trend = FALSE,
        spatial_only = FALSE,
        control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
      )
    
    mi_model <-
      sdmTMB(
        formula = mi ~ -1 + log_depth_scaled + log_depth_scaled2
        +  as.factor(year),
        data = dat,
        time = "year",
        spde = c_spde,
        anisotropy = TRUE,
        silent = TRUE,
        spatial_trend = FALSE,
        spatial_only = FALSE,
        control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
      )
    
    po2_model <-
      sdmTMB(
        formula = po2 ~ -1 + log_depth_scaled + log_depth_scaled2+as.factor(year),
        data = dat,
        time = "year",
        spde = c_spde,
        anisotropy = TRUE,
        silent = TRUE,
        spatial_trend = FALSE,
        spatial_only = FALSE,
        control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
      )
    
    o2_model <-
      sdmTMB(
        formula = o2 ~ -1 + log_depth_scaled + log_depth_scaled2  +  as.factor(year),
        data = dat,
        time = "year",
        spde = c_spde,
        anisotropy = TRUE,
        silent = TRUE,
        spatial_trend = FALSE,
        spatial_only = FALSE,
        control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
      )
    
    # Now make predictions onto grid
    wcvi_grid$loc = seq(1, nrow(wcvi_grid))
    wcvi_grid$log_depth_scaled <- (log(wcvi_grid$depth) - mean(log(dat$depth))) / sd(log(dat$depth))
    wcvi_grid$log_depth_scaled2 <- wcvi_grid$log_depth_scaled^2
    # truncate limits based on haul filters for OR above
    df = expand.grid(loc = unique(wcvi_grid$loc),
                     year = unique(dat$year))
    
    df = left_join(df, wcvi_grid)
    df$log_depth_scaled = matrix(df$log_depth_scaled, ncol = 1)
    df$log_depth_scaled2 = matrix(df$log_depth_scaled2, ncol = 1)
    
    pred_temp = predict(temp_model,
                        newdata = df,
                        return_tmb_object = F)
    pred_temp$unscaled_est <-
      pred_temp$est * sd(dat$temp_unscaled, na.rm = T) + mean(dat$temp_unscaled, na.rm = T)
    
    pred_mi <- predict(mi_model,
                       newdata = df,
                       return_tmb_object = F)
    pred_mi$unscaled_est <-
      pred_mi$est * sd(dat$mi_unscaled, na.rm = T) + mean(dat$mi_unscaled, na.rm = T)
    pred_po2 <- predict(po2_model,
                        newdata = df,
                        return_tmb_object = F)
    pred_po2$unscaled_est <-
      pred_po2$est * sd(dat$po2_unscaled, na.rm = T) + mean(dat$po2_unscaled, na.rm = T)
    
    pred_o2 <- predict(o2_model,
                       newdata = df,
                       return_tmb_object = F)
    pred_o2$unscaled_est <-
      pred_o2$est * sd(dat$o2_unscaled, na.rm = T) + mean(dat$o2_unscaled, na.rm = T)
    
    saveRDS(pred_temp, file = "output/bc/temp.rds")
    saveRDS(pred_po2, file = "output/bc/po2.rds")
    saveRDS(pred_mi, file = "output/bc/mi.rds")
    saveRDS(pred_o2, file = "output/bc/o2.rds")
    
    saveRDS(temp_model, file = "output/bc/temp_model.rds")
    saveRDS(po2_model, file = "output/bc/po2_model.rds")
    saveRDS(o2_model, file = "output/bc/o2_model.rds")
    saveRDS(mi_model, file = "output/bc/mi_model.rds")
  }
  
  if (!fit_env) {
    # Now predict sablefish catch rates
    # load sablefish catch model (best performing model)
    # Now make predictions onto grid
    wcvi_grid$loc = seq(1, nrow(wcvi_grid))
    wcvi_grid$log_depth_scaled <- (log(wcvi_grid$depth) - mean(log(dat$depth))) / sd(log(dat$depth))
    wcvi_grid$log_depth_scaled2 <- wcvi_grid$log_depth_scaled^2
    # truncate limits based on haul filters for OR above
    df = expand.grid(loc = unique(wcvi_grid$loc),
                     year = unique(dat$year))
    
    df = left_join(df, wcvi_grid)
    df$log_depth_scaled = matrix(df$log_depth_scaled, ncol = 1)
    df$log_depth_scaled2 = matrix(df$log_depth_scaled2, ncol = 1)
    sablefish_model <- readRDS("output/bc/model_9_MI.rds")
    
    pred_temp <- readRDS("output/bc/temp.rds")
    pred_mi <- readRDS("output/bc/mi.rds")
    pred_po2 <- readRDS("output/bc/po2.rds")
    pred_o2 <- readRDS("output/bc/o2.rds")
   
    # Add predicted (scaled) temp, mi and po2 to grid
    # Now make predictions onto grid
    
    
    df$temp <- matrix(pred_temp$est)
    df$po2 <- matrix(pred_po2$est)
    df$mi <- matrix(pred_mi$est)
    df$o2 <- matrix(pred_o2$est)
    
    sablefish_pred <-
      predict(sablefish_model,
              newdata = df,
              return_tmb_object = F)
    
 
    saveRDS(sablefish_pred, file = "output/bc/catch_rate.rds")
  }
}


if (!run_models) {
  pred_temp <- readRDS("output/bc/temp.rds")
  pred_mi <- readRDS("output/bc/mi.rds")
  pred_po2 <- readRDS("output/bc/po2.rds")
  pred_o2 <- readRDS("output/bc/o2.rds")
  sablefish_pred_po2 <- readRDS("output/bc/catch_rate.rds")
}
pred_all <- dplyr::filter(pred_temp, year>=2010)
pred_po2 <- dplyr::filter(pred_po2, year>=2010)
pred_o2 <- dplyr::filter(pred_o2, year>=2010)
pred_mi <- dplyr::filter(pred_mi, year>=2010)
pred_all$temperature <- pred_all$est_rf
pred_all$mi <- pred_mi$est_rf
pred_all$po2 <- pred_po2$est_rf
pred_all$o2 <- pred_o2$est_rf

# Plot output
# A short function for plotting our predictions:
plot_map <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap( ~ year) +
    coord_fixed()
}

plot_map2 <- function(dat, column = "prediction") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap( ~ year, ncol = 6) +
    coord_fixed()
}

plot_map2(pred_all, column = "o2")
