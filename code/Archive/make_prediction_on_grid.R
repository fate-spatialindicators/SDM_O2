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

run_models <- T
fit_env <- T

# load models and data
formulas <- get_models()
dat <- load_data()

if (fit_env)
  dat <- load_all_hauls()
if (!fit_env)
  dat <- load_data()


wc_grid <- readRDS("survey_data/wc_grid.rds")
dat$log_depth_scaled <- scale(log(dat$depth))
dat$log_depth_scaled2 <- dat$log_depth_scaled ^ 2
dat$jday_scaled <- scale(dat$jday)
dat$jday_scaled2 <- dat$jday_scaled ^ 2
dat$mi_unscaled <- as.numeric(dat$mi)
dat$po2_unscaled <- as.numeric(dat$po2)
dat$o2_unscaled <- as.numeric(dat$o2)
dat$temp_unscaled <- as.numeric(dat$temp)

dat$mi <- as.numeric(scale(dat$mi))
dat$o2 <- as.numeric(scale(dat$o2))
dat$po2 <- as.numeric(scale(dat$po2))
dat$temp <- as.numeric(scale(dat$temp))

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
        +  as.factor(year) + jday_scaled + jday_scaled2,
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
        +  as.factor(year) + jday_scaled + jday_scaled2,
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
        formula = po2 ~ -1 + log_depth_scaled + log_depth_scaled2+as.factor(year) + jday_scaled + jday_scaled2,
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
        formula = o2 ~ -1 + log_depth_scaled + log_depth_scaled2  +  as.factor(year) + jday_scaled + jday_scaled2,
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
    wc_grid$loc = seq(1, nrow(wc_grid))
    
    # truncate limits based on haul filters for OR above
    df = expand.grid(loc = unique(wc_grid$loc),
                     year = unique(dat$year))
    df = left_join(df, wc_grid)
    df$jday_scaled = 0
    df$jday_scaled2 = 0
    df$log_depth_scaled = matrix(df$log_depth_scaled, ncol = 1)
    df$log_depth_scaled2 = matrix(df$log_depth_scaled2, ncol = 1)
    df$jday_scaled = matrix(df$jday_scaled, ncol = 1)
    df$jday_scaled2 = matrix(df$jday_scaled2, ncol = 1)
    
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
    
    saveRDS(pred_temp, file = "output/wc/temp.rds")
    saveRDS(pred_po2, file = "output/wc/po2.rds")
    saveRDS(pred_mi, file = "output/wc/mi.rds")
    saveRDS(pred_o2, file = "output/wc/o2.rds")
    
    saveRDS(temp_model, file = "output/wc/temp_model.rds")
    saveRDS(po2_model, file = "output/wc/po2_model.rds")
    saveRDS(o2_model, file = "output/wc/o2_model.rds")
    saveRDS(mi_model, file = "output/wc/mi_model.rds")
  }
  
  if (!fit_env) {
    # Now predict sablefish catch rates
    # load sablefish catch model (best performing model)
    wc_grid <- readRDS("survey_data/wc_grid.rds")
    sablefish_model_po2 <- readRDS("output/wc/model_13_MI.rds")
    
    pred_temp <- readRDS("output/wc/temp.rds")
    pred_mi <- readRDS("output/wc/mi.rds")
    pred_po2 <- readRDS("output/wc/po2.rds")
    pred_o2 <- readRDS("output/wc/o2.rds")
    pred_temp <- dplyr::filter(pred_temp, year >= 2010)
    pred_mi <- dplyr::filter(pred_mi, year >= 2010)
    pred_po2 <- dplyr::filter(pred_po2, year >= 2010)
    pred_o2 <- dplyr::filter(pred_o2, year >= 2010)
    # Add predicted (scaled) temp, mi and po2 to grid
    # Now make predictions onto grid
    wc_grid$loc = seq(1, nrow(wc_grid))
    
    # truncate limits based on haul filters for OR above
    df = expand.grid(loc = unique(wc_grid$loc),
                     year = unique(dat$year))
    df <- dplyr::filter(df, year >=2010)
    
    df = left_join(df, wc_grid)
    df$jday_scaled = 0
    df$jday_scaled2 = 0
    df$log_depth_scaled = matrix(df$log_depth_scaled, ncol = 1)
    df$log_depth_scaled2 = matrix(df$log_depth_scaled2, ncol = 1)
    df$jday_scaled = matrix(df$jday_scaled, ncol = 1)
    df$jday_scaled2 = matrix(df$jday_scaled2, ncol = 1)
    
    df$temp <- pred_temp$est
    df$po2 <- pred_po2$est
    df$mi <- pred_mi$est
    df$o2 <- pred_o2$est
    
    sablefish_pred_po2 <-
      predict(sablefish_model_po2,
              newdata = df,
              return_tmb_object = F)
    
 
    saveRDS(sablefish_pred_po2, file = "output/wc/catch_rate.rds")
  }
}


if (!run_models) {
  pred_temp <- readRDS("output/wc/temp.rds")
  pred_mi <- readRDS("output/wc/mi.rds")
  pred_po2 <- readRDS("output/wc/po2.rds")
  pred_o2 <- readRDS("output/wc/o2.rds")
  sablefish_pred_po2 <- readRDS("output/wc/catch_rate.rds")
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

p_temp <- plot_map2(pred_all, "temperature") +
  scale_fill_viridis_c(option = "viridis")
p_o2 <- plot_map2(pred_all, "o2") +
  scale_fill_viridis_c(option = "viridis")
p_po2 <- plot_map2(pred_all, "po2") +
  scale_fill_viridis_c(option = "viridis")
p_mi <- plot_map2(pred_all, "mi") +
  scale_fill_viridis_c(option = "viridis")

grid.arrange(p_temp, p_o2, nrow = 2, ncol = 1)
grid.arrange(p_po2, p_mi, nrow = 2, ncol = 1)

plot_map2(pred, "prediction") +
  scale_fill_viridis_c(option = "viridis")



plot_map(pred_temp, "unscaled_est") +
  scale_fill_viridis_c(option = "viridis") +
  ggtitle("Predicted Temp") +
  labs("T")

plot_map(pred_po2, "unscaled_est") +
  scale_fill_viridis_c(option = "viridis") +
  ggtitle("Predicted po2") +
  labs("pO2")

plot_map(pred_mi, "unscaled_est") +
  scale_fill_viridis_c(option = "viridis") +
  ggtitle("Predicted MI") +
  labs("MI")

# Make a map of the difference in estimated SCALED po2 and MI
delta_mi_po2 <- pred_mi
delta_mi_po2$deltao2 <- pred_mi$est - pred_po2$est

plot_map(delta_mi_po2, "deltao2") +
  scale_fill_viridis_c(option = "plasma") +
  ggtitle("Scaled MI - Scaled po2") +
  labs("difference")


plot_map(sablefish_pred_po2, "exp(est)") +
  scale_fill_viridis_c(option = "viridis") +
  ggtitle("Predicted catch rate") +
  labs("Catch Rate (kg / km2)")




## Fit other models and compare predictions spatially

wc_grid$loc = seq(1, nrow(wc_grid))

# truncate limits based on haul filters for OR above
df = expand.grid(loc = unique(wc_grid$loc),
                 year = unique(dat$year))
df = left_join(df, wc_grid)
df$jday_scaled = 0
df$jday_scaled2 = 0
df$log_depth_scaled = matrix(df$log_depth_scaled, ncol = 1)
df$log_depth_scaled2 = matrix(df$log_depth_scaled2, ncol = 1)
df$jday_scaled = matrix(df$jday_scaled, ncol = 1)
df$jday_scaled2 = matrix(df$jday_scaled2, ncol = 1)
df$temp <- pred_temp$est
df$po2 <- pred_po2$est
df$mi <- pred_mi$est

formula <- paste0("cpue_kg_km2 ~ -1 +", formulas[14])
sablefish_model_mi <- sdmTMB(
  formula = as.formula(formula),
  data = dat,
  time = NULL,
  reml = TRUE,
  spde = c_spde,
  family = tweedie(link = "log"),
  anisotropy = TRUE,
  spatial_only = TRUE,
  silent = TRUE
)
sablefish_pred_mi <-
  predict(sablefish_model_mi,
          newdata = df,
          return_tmb_object = F)

po2_vs_mi <- sablefish_pred_mi
po2_vs_mi$delta_yhat <-
  exp(sablefish_pred_mi$est) - exp(sablefish_pred_po2$est)
plot_map(po2_vs_mi, "delta_yhat") +
  scale_fill_viridis_c(option = "viridis") +
  ggtitle("Difference in predicted catch rate (mi - po2 model)") +
  labs("Catch Rate (kg / km2)")

# make another prediction assuming all po2 is above breakpoint
# need to load the model
sablefish_model_po2 <- readRDS("output/wc/model_13_MI.rds")

bp_pars <- get_bp_parameters(sablefish_model_po2)
df2 <- df

df2$po2 <- rep(bp_pars[1], times = nrow(df2))

sablefish_pred_highpo2 <-
  predict(sablefish_model_po2,
          newdata = df2,
          return_tmb_object = F)

delta_predict <- sablefish_pred_po2
delta_predict$est <-
  sablefish_pred_po2$est - sablefish_pred_highpo2$est
plot_map(delta_predict, "exp(est)") +
  scale_fill_viridis_c(option = "viridis") +
  ggtitle("Predicted po2 effect")


# Plot po2 anomolies
pred_ave <- pred_po2 %>%
  group_by(loc) %>%
  summarise(Ave_po2 = mean(est_rf))

# big ass loop
pred_po2$anomoly <- rep(NA, nrow(pred_po2))
for (i in 1:nrow(pred_ave)) {
  index <- which(pred_po2$loc == pred_ave$loc[i])
  anomoly <- pred_po2$est_rf[index] - pred_ave$Ave_po2[i]
  pred_po2$anomoly[index] <- anomoly
}

plot_map2(pred_po2, "est_rf") +
  scale_fill_viridis_c(option = "viridis")
