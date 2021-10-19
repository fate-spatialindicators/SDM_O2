
source("code/mi_functions.R")
library(raster)
library(dplyr)
library(sp)
library(gsw)
library(rgdal)
library(future)
library(ggplot2)
library(viridis)

spc <- "sablefish"
fit.model <- F # do you want to fit SDM to environmental variables and impute missing values?
years <- 2010:2015 # designate years to use, must be 2010:2015 if using trawl-based data
constrain_latitude <- F # do you want to constraint trawl data N of 43 degrees latitude
compare_sources <- F # do you want to run models comparing trawl and J-SCOPE covariates?
no_depth <- FALSE # Do you want to run models w/out a depth effect?
use_cv = FALSE # specify whether to do cross validation or not
use_AIC = TRUE # specify whether to use AIC
use_jscope <- F # specify whether to only use J-SCOPE based estimates.  Overrides compare_sources and fit.model
dat <- load_data(fit.model, spc, constrain_latitude)
# only look at intermediate sized sablefish, 0.5 - 6 kg)
dat$cpue_kg_km2 <- dat$cpue_kg_km2 * (dat$p2 + dat$p3) 


kelvin = 273.15
kb <- boltz <- 0.000086173324
Eo <- 0.4525966 
Ao <- 5.780791e-06
n <- -0.303949 
B = 1000 # size in grams, roughly average
avebn <- 0.1124206
logAo <- log(Ao)

t.range <- 2:12
k.range <- t.range + 273.15

plot.phi <- data.frame(inv_temp =1/(k.range * kb), 
                       logbnpo2 = -logAo - Eo/(kb * k.range))



plot.obs <- data.frame(x = 1/((kelvin + 12) * kb), y = - logAo - Eo/(kb* (kelvin + 12)))
dat$inv_temp <- 1 / (boltz * (dat$temp+kelvin))
dat$logbnpo2 <- log(avebn * dat$po2)

# only look at depths that are commonly inhabited by sablefish
dat <- dplyr::filter(dat, depth >=100)

pointsize <- 0.5
alpha2use <- 0.7
zlims <- c(4,8)

po2plot <- ggplot() + 
  geom_point(data = dat, aes(x = po2, y = -(depth), col = log(cpue_kg_km2)), alpha = alpha2use,size = pointsize) +
  scale_x_continuous() +
  scale_colour_viridis(limits = zlims,oob = scales::squish,name = bquote('log(biomass)'~(kg~km^2))) +
  labs(x = bquote(pO[2]~"(atm)"), y = "Depth (m)") +
  #geom_line(data = plot.phi, aes(x = inv_temp, y = logbnpo2), size = 1.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(legend.position = c(0.625, 0.35),
            legend.direction = "vertical")+
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text = element_text(size = 12))
  

miplot <- ggplot() + 
  geom_point(data = dat, aes(x = mi, y = -(depth), col = log(cpue_kg_km2)), alpha = alpha2use, size = pointsize, show.legend = FALSE) +
  scale_x_continuous() +
  scale_colour_viridis(limits = zlims,oob = scales::squish,name = bquote('log(biomass)'~(kg~km^2))) +
  labs(x = "Metabolic Index", y = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text = element_text(size = 12))

templot <- ggplot() + 
  geom_point(data = dat, aes(x = temp, y = -(depth), col = log(cpue_kg_km2)), alpha = alpha2use, size = pointsize,show.legend= FALSE) +
  scale_x_continuous() +
  scale_colour_viridis(limits = zlims,oob = scales::squish, name = bquote('log(biomass)'~(kg~km^2))) +
  labs(x = "Temperature", y = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text = element_text(size = 12))


multipanel <- gridExtra::grid.arrange(templot, po2plot, miplot, nrow = 1, ncol = 3)

ggsave("plots/env_vscatch.png", multipanel, width = 8, height = 3.5, units = "in", device = "png")




