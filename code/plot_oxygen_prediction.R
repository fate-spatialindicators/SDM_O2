# load handy functions
source("code/mi_functions.R")

# Need to load the data to get scaled stuff
fit.model <- F # do you want to fit SDM to environmental variables and impute missing values?
years <- 2010:2015 # designate years to use, must be 2010:2015 if using trawl-based data
constrain_latitude <- F # do you want to constraint trawl data N of 43 degrees latitude
compare_sources <- F # do you want to run models comparing trawl and J-SCOPE covariates?
no_depth <- FALSE # Do you want to run models w/out a depth effect?
use_cv = FALSE # specify whether to do cross validation or not
use_AIC = TRUE # specify whether to use AIC
use_jscope <- F # specify whether to only use J-SCOPE based estimates.  Overrides compare_sources and fit.model

dat <- load_data(spc = "sablefish", constrain_latitude, fit.model)
# rescale variables
mean.depth <- mean(log(dat$depth))
sd.depth <- sd(log(dat$depth))

dat$log_depth_scaled <- as.numeric(scale(log(dat$depth)))
dat$log_depth_scaled2 <- dat$log_depth_scaled^2
dat$jday_scaled <- as.numeric(scale(dat$julian_day))
dat$jday_scaled2 <- as.numeric(scale(log(dat$julian_day) ^ 2))

mean.po2 <- mean(dat$po2)
std.po2 <- sd(dat$po2)
mean.do <- mean(dat$o2)
std.do <- sd(dat$o2)
mean.mi <- mean(dat$mi)
std.mi <- sd(dat$mi)

# use this to fit a spatio-temporal model of pO2, so must scale it first
dat$po2 <- as.numeric(scale(dat$po2))

dat$X <- dat$longitude
dat$Y <- dat$latitude

# Load best modelfor sablefish
m <- readRDS("output/wc/model_8_sablefish.rds")
c_spde <-make_mesh(data = dat, xy_cols = c("X", "Y"), n_knots = 250) # choose # knots


# removed jday_scaled2 because may not have converged
po2_model <-  sdmTMB(formula = po2 ~ -1 + log_depth_scaled + log_depth_scaled2 
                     +  as.factor(year) + jday_scaled + jday_scaled2,
                     data = dat,
                     time = "year", spde = c_spde, anisotropy = TRUE,
                     silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                     control = sdmTMBcontrol(step.min = 0.01, step.max = 1))

wc_grid$loc = seq(1,nrow(wc_grid))

# truncate limits based on haul filters for OR above
df = expand.grid(loc=unique(wc_grid$loc),
                 year = unique(dat$year))
df = left_join(df,wc_grid)
df$jday_scaled = 0
df$jday_scaled2 = 0
df$log_depth_scaled = matrix((log(-df$depth) - mean.depth)/sd.depth, ncol=1)
df$log_depth_scaled2 = matrix(df$log_depth_scaled^2, ncol=1)

pred_po2 <- predict(po2_model,
                    newdata = df,
                    return_tmb_object = F)


m_po2 <- readRDS("output/wc/model_13_sablefish.rds")
m_o2 <- readRDS("output/wc/model_11_sablefish.rds")
m_null <- readRDS("output/wc/model_2_MI.rds")

# predict_po2 <- predict(m_po2)
# predict_o2 <- predict(m_o2)
# predict_null <- predict(m_null)
# 
# # look at difference in predictions between oxygen models and null
# 
# predict_po2$delta <- predict_po2$est - predict_null$est
# predict_o2$delta <- predict_o2$est - predict_null$est
# 
# 
# library(scales)
# ggplot(predict_po2, aes(longitude, latitude, col = delta)) + scale_colour_gradient2() +
#   geom_point(alpha=0.6) + facet_wrap(~year) + labs(col = "Delta (po2 - null)")
# ggplot(predict_o2, aes(longitude, latitude, col = delta)) + scale_colour_gradient2() +
#   geom_point(alpha=0.6) + facet_wrap(~year) + labs(col = "Delta (o2 - null)")
# 
# # just show the values below the breakpoint
# fixed_effects_po2 <- m_po2$model$par
# bp_ind <- grep("b_threshold", names(fixed_effects_po2))
# bp_po2 <- fixed_effects_po2[bp_ind[2]]
# slope_po2 <- fixed_effects_po2[bp_ind[1]]
# 
# # create new column in data that contains breakpoint if po2 is above breakpoint, and lists the actual po2 otherwise
# zeros <- rep(bp_po2, times = nrow(dat))
# dat$po2_bp <- zeros
# 
# dat$po2_bp[dat$po2<=bp_po2] = back.convert(dat$po2[dat$po2<=bp_po2], mean.po2, std.po2)
# 
# # plot the po2 only below the threshold
# ggplot(dat, aes(longitude, latitude, col = po2_bp)) + scale_colour_gradient2() +
#   geom_point(alpha=0.6) + facet_wrap(~year) + labs(col = "p02 below breakpoint")
# 
# # plot the effect size from po2
# above <- rep(exp(bp_po2 * slope_po2), times = nrow(dat))
# dat$po2_effect <- above
# 
# dat$po2_effect[dat$po2<=bp_po2] = exp((dat$po2[dat$po2<=bp_po2])* slope_po2)
# ggplot(dat, aes(longitude, latitude, col = po2_effect)) + scale_colour_gradient2() +
#   geom_point(alpha=0.6) + facet_wrap(~year)
# 
# # plot the log(cpue)
# 
# ggplot(dat, aes(longitude, latitude, col = log(cpue_kg_km2 ))) + scale_colour_gradient2() +
#   geom_point(alpha=0.6) + facet_wrap(~year)
