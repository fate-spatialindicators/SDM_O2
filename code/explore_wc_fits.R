# partial residual plots of MI, o2 and po2

# Load functions and packages
source("code/mi_functions.R")
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)

# load models
formulas <- get_models()
dat <- load_data()


par(mfrow=c(1,3))
model_numbers <- c(10, 12, 14)
predictors <- c("o2", "po2", "mi")

for (j in 1:length(model_numbers)) {
model_number <- model_numbers[j]
model_name <- paste0("output/wc/model_",model_number,"_MI.rds")
m <- readRDS(file = model_name)
x <- dat[,which(colnames(dat)==predictors[j])]
partial_residual<- get_partial_resid_bp(m =m , x = x)
plot(x = x, y = partial_residual,
     type = "p",
     pch = 21,
     bg = "black",
     xlab = predictors[j],
     ylab = "Partial Residual"
)
}


par(mfrow = c(1,1))
model_number <-10
model_name <- paste0("output/wc/model_",model_number,"_MI.rds")
m <- readRDS(file = model_name)

# compare fitted depth effects under the o2 versus base model. 
# depth coefficients
depth_coef <- m$sd_report$par.random[1:2]

depth <- dat$depth
s_depth <- scale(depth)

min.depth <- min(s_depth)
max.depth <- max(s_depth)
depth.list <- seq(from = min(s_depth), to = max(s_depth), length.out = 100)

depth.effect <-rep(NA, times = length(depth.list))

for (i in 1:length(depth.list))  depth.effect[i] <- depth_coef[1] * depth.list[i] + depth_coef[2] * I(depth.list[i]^2)

plot(depth.list, depth.effect, 
     type = "l", 
     lwd = 2,
     xlab = "depth",
     ylab = "effect size")

# repeat with null model

model_number <- 1
model_name <- paste0("output/wc/model_",model_number,"_MI.rds")
m <- readRDS(file = model_name)

# compare fitted depth effects under the o2 versus base model. 
# depth coefficients
depth_coef <- m$sd_report$par.random[1:2]

depth <- dat$depth
s_depth <- scale(depth)

min.depth <- min(s_depth)
max.depth <- max(s_depth)
depth.list <- seq(from = min(s_depth), to = max(s_depth), length.out = 100)

depth.effect <-rep(NA, times = length(depth.list))

for (i in 1:length(depth.list))  depth.effect[i] <- depth_coef[1] * depth.list[i] + depth_coef[2] * I(depth.list[i]^2)

lines(depth.list, depth.effect,
      type = "l",
      lwd = 2,
      col = "darkblue")
