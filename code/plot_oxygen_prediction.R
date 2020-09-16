# load handy functions
source("code/mi_functions.R")

# Script to compare fitted Oxygen models
dat <- load_data()
# rescale variables
dat$depth = scale(log(dat$depth))
dat$temp = scale(dat$temp)

mean.po2 <- mean(dat$po2)
std.po2 <- sd(dat$po2)
mean.do <- mean(dat$o2)
std.do <- sd(dat$o2)
mean.mi <- mean(dat$mi)
std.mi <- sd(dat$mi)
dat$mi = scale(dat$mi)
dat$o2 = scale(dat$o2)
dat$po2 <- scale(dat$po2)


# first load models
m_po2 <- readRDS("output/wc/model_13_MI.rds")
m_o2 <- readRDS("output/wc/model_11_MI.rds")
m_null <- readRDS("output/wc/model_2_MI.rds")

predict_po2 <- predict(m_po2)
predict_o2 <- predict(m_o2)
predict_null <- predict(m_null)

# look at difference in predictions between oxygen models and null

predict_po2$delta <- predict_po2$est - predict_null$est
predict_o2$delta <- predict_o2$est - predict_null$est


library(scales)
ggplot(predict_po2, aes(longitude, latitude, col = delta)) + scale_colour_gradient2() +
  geom_point(alpha=0.6) + facet_wrap(~year) + labs(col = "Delta (po2 - null)")
ggplot(predict_o2, aes(longitude, latitude, col = delta)) + scale_colour_gradient2() +
  geom_point(alpha=0.6) + facet_wrap(~year) + labs(col = "Delta (o2 - null)")

# just show the values below the breakpoint
fixed_effects_po2 <- m_po2$model$par
bp_ind <- grep("b_threshold", names(fixed_effects_po2))
bp_po2 <- fixed_effects_po2[bp_ind[2]]
slope_po2 <- fixed_effects_po2[bp_ind[1]]

# create new column in data that contains breakpoint if po2 is above breakpoint, and lists the actual po2 otherwise
zeros <- rep(bp_po2, times = nrow(dat))
dat$po2_bp <- zeros

dat$po2_bp[dat$po2<=bp_po2] = back.convert(dat$po2[dat$po2<=bp_po2], mean.po2, std.po2)

# plot the po2 only below the threshold
ggplot(dat, aes(longitude, latitude, col = po2_bp)) + scale_colour_gradient2() +
  geom_point(alpha=0.6) + facet_wrap(~year) + labs(col = "p02 below breakpoint")

# plot the effect size from po2
above <- rep(exp(bp_po2 * slope_po2), times = nrow(dat))
dat$po2_effect <- above

dat$po2_effect[dat$po2<=bp_po2] = exp((dat$po2[dat$po2<=bp_po2])* slope_po2)
ggplot(dat, aes(longitude, latitude, col = po2_effect)) + scale_colour_gradient2() +
  geom_point(alpha=0.6) + facet_wrap(~year)

# plot the log(cpue)

ggplot(dat, aes(longitude, latitude, col = log(cpue_kg_km2 ))) + scale_colour_gradient2() +
  geom_point(alpha=0.6) + facet_wrap(~year)
