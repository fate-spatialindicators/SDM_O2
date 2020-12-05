# partial residual plots of MI, o2 and po2

# Load functions and packages
source("code/mi_functions.R")
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)
library(ggplot2)

# load models
formulas <- get_models()
dat <- load_data()
model_numbers <- c(10, 12, 14)
predictors <- c("o2", "po2", "mi")
model_number <- 12

#for (j in 1:length(model_numbers)) {
#x <- dat[,which(colnames(dat)==predictors[2])]

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
dat$temp <- as.numeric(scale(dat$po2))

c_spde <-make_mesh(data = dat, xy_cols = c("X", "Y"), n_knots = 250) # choose # knots


# Make grid of environmental conditions
temp_model <- sdmTMB(formula = temp ~ -1 + log_depth_scaled + log_depth_scaled2 
                     +  as.factor(year) + jday_scaled + jday_scaled2,
                     data = dat,
                     time = "year", spde = c_spde, anisotropy = TRUE,
                     silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                     control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

mi_model <-  sdmTMB(formula = mi ~ -1 + log_depth_scaled + log_depth_scaled2 
                    +  as.factor(year) + jday_scaled + jday_scaled2,
                    data = dat,
                    time = "year", spde = c_spde, anisotropy = TRUE,
                    silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                    control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

# removed jday_scaled2 because may not have converged
po2_model <-  sdmTMB(formula = po2 ~ -1 + log_depth_scaled + log_depth_scaled2 
                     +  as.factor(year) + jday_scaled + jday_scaled2,
                     data = dat,
                     time = "year", spde = c_spde, anisotropy = TRUE,
                     silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                     control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

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


pred_mi <- predict(mi_model,
                   newdata = df,
                   return_tmb_object = F)

pred_po2 <- predict(po2_model,
                    newdata = df,
                    return_tmb_object = F)


# Add predicted temp, mi and po2 to grid
df$temp <- matrix(pred_temp$est, ncol = 1)
df$po2 <- matrix(pred_po2$est, ncol = 1)
df$mi <- matrix(pred_mi$est, ncol = 1)

# Now predict sablefish catch rates
po2_predict <- predict(po2_model, newdata =df, return_tmb_object = F)

# Plot output
# A short function for plotting our predictions:
plot_map <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

plot_map(po2_predict, "est") +
  scale_fill_viridis_c() +
  ggtitle("Predicted log catch rate")

