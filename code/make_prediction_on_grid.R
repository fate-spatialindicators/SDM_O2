# Load functions and packages
#devtools::install_github("pbs-assess/sdmTMB")
source("code/mi_functions.R")
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)
library(ggplot2)

# load models
formulas <- get_models()
dat <- load_data()
wc_grid <- readRDS("survey_data/wc_grid.rds")

dat <- load_data()
# transform
dat$log_depth_scaled <- scale(log(dat$depth))
dat$log_depth_scaled2 <- dat$log_depth_scaled^2
dat$jday_scaled <- scale(dat$julian_day)
dat$jday_scaled2 <- dat$jday_scaled^2
dat$X <- dat$longitude
dat$Y <- dat$latitude
dat$mi_unscaled <- as.numeric(dat$mi)
dat$po2_unscaled <- as.numeric(dat$po2)
dat$o2_unscaled <- as.numeric(dat$o2)
dat$temp_unscaled <- as.numeric(dat$temp)

dat$mi <- as.numeric(scale(dat$mi))
dat$o2 <- as.numeric(scale(dat$o2))
dat$po2 <- as.numeric(scale(dat$po2))
dat$temp <- as.numeric(scale(dat$temp))

c_spde <-make_mesh(data = dat, xy_cols = c("X", "Y"), n_knots = 250) # choose # knots


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
    formula = po2 ~ -1 + log_depth_scaled + log_depth_scaled2  +  as.factor(year) + jday_scaled + jday_scaled2,
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
wc_grid$loc = seq(1,nrow(wc_grid))

# truncate limits based on haul filters for OR above
df = expand.grid(loc=unique(wc_grid$loc),
                 year = unique(dat$year))
df = left_join(df,wc_grid)
df$jday_scaled = 0
df$jday_scaled2 = 0
df$log_depth_scaled = matrix(df$log_depth_scaled, ncol=1)
df$log_depth_scaled2 = matrix(df$log_depth_scaled2, ncol=1)
df$jday_scaled = matrix(df$jday_scaled, ncol=1)
df$jday_scaled2 = matrix(df$jday_scaled2, ncol=1)

pred_temp = predict(temp_model,
                    newdata=df,
                    return_tmb_object = F)
pred_temp$unscaled_est <- pred_temp$est *sd(dat$temp_unscaled, na.rm = T) + mean(dat$temp_unscaled, na.rm = T) 

pred_mi <- predict(mi_model,
                   newdata = df,
                   return_tmb_object = F)
pred_mi$unscaled_est <- pred_mi$est *sd(dat$mi_unscaled, na.rm = T) + mean(dat$mi_unscaled, na.rm = T) 
pred_po2 <- predict(po2_model,
                    newdata = df,
                    return_tmb_object = F)
pred_po2$unscaled_est <- pred_po2$est *sd(dat$po2_unscaled, na.rm = T) + mean(dat$po2_unscaled, na.rm = T) 

# Add predicted (scaled) temp, mi and po2 to grid
df$temp <- pred_temp$est
df$po2 <- pred_po2$est
df$mi <- pred_mi$est

# Now predict sablefish catch rates
# Fit sablefish catch model (best performing model)
formula <- paste0("cpue_kg_km2 ~ -1 +", formulas[13])
sablefish_model <-sdmTMB(formula = as.formula(formula),
                         data = dat,
                         time = NULL,
                         reml = TRUE,
                         spde = c_spde,
                         family = tweedie(link = "log"),
                         anisotropy = TRUE,
                         spatial_only =TRUE,
                         silent = TRUE)


sablefish_predict <- predict(sablefish_model, newdata =df, return_tmb_object = F)

# Plot output
# A short function for plotting our predictions:
plot_map <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}


plot_map(pred_po2, "unscaled_est") +
scale_fill_viridis_c(option = "plasma") +
  ggtitle("Predicted po2") +
  labs("pO2")

plot_map(pred_temp, "unscaled_est") +
  scale_fill_viridis_c(option = "plasma") +
  ggtitle("Predicted Temp") +
  labs("T")


plot_map(pred_mi, "unscaled_est") +
  scale_fill_viridis_c(option = "plasma") +
  ggtitle("Predicted MI") +
  labs("MI")


plot_map(sablefish_predict, "exp(est)") +
  scale_fill_viridis_c(option = "viridis") +
  ggtitle("Predicted catch rate") +
 labs("Catch Rate (kg / km2)")

# make another prediction assuming all po2 is above breakpoint
bp_pars <- get_bp_parameters(sablefish_model)
df2 <- df

df2$po2 <- rep(bp_pars[1], nrow(df2))
sablefish_predict_highpo2 <- predict(sablefish_model, newdata =df2, return_tmb_object = F)

delta_predict <- sablefish_predict
delta_predict$est <- sablefish_predict$est - sablefish_predict_highpo2$est
plot_map(delta_predict, "exp(est)") +
  scale_fill_viridis_c(option = "viridis") +
  ggtitle("Predicted po2 effect")
