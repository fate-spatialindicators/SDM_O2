library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)
library(tweedie)

# load handy functions
source("code/mi_functions.R")
m_po2 <- readRDS("output/wc/model_7_sablefish.rds") 

phi <- exp(m_po2$model$par[6])
theta <- m_po2$model$par[5]
power <- 1 + exp(theta) / (1 + exp(theta))


dat <- load_data(fit.model=F, spc="sablefish", constrain_latitude= F)
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
dat$year <- as.factor(dat$year)
ln_yhat <- predict(m_po2, dat)
dat$est <- exp(ln_yhat$est)


n.data <- nrow(dat)
n.sims <- 1000
sim.output <- matrix(NA, nrow = n.data, ncol = n.sims)

for (i in 1:n.sims) sim.output[,i] <- rtweedie(n.data, mu = dat$est, phi = phi, power = power)
observed_quantiles <- rep(NA, times = n.data)
for (i in 1:n.data) observed_quantiles[i] <- length(which(sim.output[i,]<=dat$cpue_kg_km2[i])) / n.sims

# make qq plot
png(filename = "plots/qqplot.png", width = 6, height = 5, units = "in", res = 150)
plot((1:n.data) / n.data, y= sort(observed_quantiles, decreasing = F),
     type = "p",
     pch = 21,
     bg = "black",
     cex = 0.75,
     xlab = "Theoretical Quantiles",
     ylab = "Observed Quantiles",
     xlim = c(0,1),
     ylim = c(0,1),
     xaxs = "i",
     yaxs = "i",
     las = 1,
     cex.lab = 1.5,
     cex.axis = 1.25)
dev.off()


