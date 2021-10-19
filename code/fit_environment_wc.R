## Code to fit spatio-temporal model of environmental characteristics 

# Install latest version of sdmTMB ----------------------------------------
#devtools::install_github("pbs-assess/sdmTMB")

# Initialize with packages and functions ----------------------------

rm(list = ls())
source("code/mi_functions.R")
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)
library(ggplot2)
library(viridis)

# load and scale data -----------------------------------------------------

wc_grid <- readRDS("survey_data/wc_grid.rds")

dat <- load_data(spc = "sablefish", constrain_latitude=F, fit.model= F)
# transform
dat$log_depth_scaled <- scale(log(dat$depth))
dat$log_depth_scaled2 <- scale(log(dat$depth) ^ 2)
dat$jday_scaled <- scale(dat$julian_day)
dat$jday_scaled2 <- scale(log(dat$julian_day) ^ 2)
dat$X <- dat$longitude
dat$Y <- dat$latitude
dat$temp <- scale(dat$temp)
dat$mi <- scale(dat$mi)
dat$o2 <- scale(dat$o2)
dat$po2 <- scale(dat$po2)

c_spde <-make_mesh(data = dat, xy_cols = c("X", "Y"), n_knots = 250) # choose # knots



temp_model <- sdmTMB(formula = temp ~ -1 + log_depth_scaled + log_depth_scaled2 
                     +  as.factor(year) + jday_scaled + jday_scaled2,
                     data = dat,
                     time = "year", spde = c_spde, anisotropy = TRUE,
                     silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                     control = sdmTMBcontrol(step.min = 0.01, step.max = 1))
temp_model_reduced <- sdmTMB(formula = temp ~ -1 + log_depth_scaled + log_depth_scaled2 
                             +  as.factor(year) + jday_scaled + jday_scaled2,
                             data = dat,
                             spde = c_spde, anisotropy = TRUE,
                             silent = TRUE, spatial_trend = FALSE, spatial_only = TRUE,
                             control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

temp_model_noyear <- sdmTMB(formula = temp ~ -1 + log_depth_scaled + log_depth_scaled2 
                              + jday_scaled + jday_scaled2,
                             data = dat,
                             spde = c_spde, anisotropy = TRUE,
                             silent = TRUE, spatial_trend = FALSE, spatial_only = TRUE,
                             control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

temp_aic <- c(AIC(temp_model_noyear), AIC(temp_model_reduced), AIC(temp_model))
temp_daic <- temp_aic - min(temp_aic)

mi_model <-  sdmTMB(formula = mi ~ -1 + log_depth_scaled + log_depth_scaled2 
                    +  as.factor(year) + jday_scaled + jday_scaled2,
                    data = dat,
                    time = "year", spde = c_spde, anisotropy = TRUE,
                    silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                    control = sdmTMBcontrol(step.min = 0.01, step.max = 1))


mi_model_reduced <-  sdmTMB(formula = mi ~ -1 + log_depth_scaled + log_depth_scaled2 
                    +  as.factor(year) + jday_scaled + jday_scaled2,
                    data = dat,
                    time = "year", spde = c_spde, anisotropy = TRUE,
                    silent = TRUE, spatial_trend = FALSE, spatial_only = TRUE,
                    control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

mi_model_noyear <-  sdmTMB(formula = mi ~ -1 + log_depth_scaled + log_depth_scaled2 
                     + jday_scaled + jday_scaled2,
                    data = dat,
                    time = "year", spde = c_spde, anisotropy = TRUE,
                    silent = TRUE, spatial_trend = FALSE, spatial_only = TRUE,
                    control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

mi_aic <- c(AIC(mi_model_noyear), AIC(mi_model_reduced), AIC(mi_model))
mi_daic <- mi_aic - min(mi_aic)

# removed jday_scaled2 because may not have converged
po2_model <-  sdmTMB(formula = po2 ~ -1 + log_depth_scaled + log_depth_scaled2 
                    +  as.factor(year) + jday_scaled + jday_scaled2,
                    data = dat,
                    time = "year", spde = c_spde, anisotropy = TRUE,
                    silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                    control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

po2_model_reduced <-  sdmTMB(formula = po2 ~ -1 + log_depth_scaled + log_depth_scaled2 
                     +  as.factor(year) + jday_scaled + jday_scaled2,
                     data = dat,
                     time = "year", spde = c_spde, anisotropy = TRUE,
                     silent = TRUE, spatial_trend = FALSE, spatial_only = TRUE,
                     control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

po2_model_noyear <-  sdmTMB(formula = po2 ~ -1 + log_depth_scaled + log_depth_scaled2 
                     + jday_scaled + jday_scaled2,
                     data = dat,
                     time = "year", spde = c_spde, anisotropy = TRUE,
                     silent = TRUE, spatial_trend = FALSE, spatial_only = TRUE,
                     control = sdmTMBcontrol(step.min = 0.01, step.max = 1))
po2_aic <- c(AIC(po2_model_noyear), AIC(po2_model_reduced), AIC(po2_model))
po2_daic <- po2_aic - min(po2_aic)

print(cbind(temp_daic, mi_daic, po2_daic))

# make prediction surfaces of environmental conditions

# temp predictions

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


pred_mi <- predict(mi_model,
                   newdata = df,
                   return_tmb_object = F)

pred_po2 <- predict(po2_model,
                    newdata = df,
                    return_tmb_object = F)
# save to output
save(file = "output/wc_temp.Rdata", pred_temp)
save(file = "output/wc_mi.Rdata", pred_mi)
save(file = "output/wc_po2.Rdata", pred_po2)


## Plot
# A short function for plotting our predictions:
plot_map <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

plot_map(pred_temp, "est") +
  scale_fill_viridis_c() +
  ggtitle("Prediction (fixed effects + all random effects)")

plot_map(pred_mi, "est") +
  scale_fill_viridis_c() +
  ggtitle("Metabolic Index")

plot_map(pred_po2, "est") +
  scale_fill_viridis_c() +
  ggtitle("pO2 ")

