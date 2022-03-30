library(TMB)
library(ggplot2)
library(viridis)
library(gsw)
library(dplyr)
library(rMR)
library(lme4)
library(RColorBrewer)
library(rstan)
library(shinystan)
library(loo)
library(glmmTMB)

library(rMR)
kb <-  8.617333262145E-5

### set up data ####
all.dat<- read.csv("data/metabolic_index_data.csv", header = T)
all.dat$po2 <- NA
#### cycle by species and measurement type, and get po2 for each ####
spc.list <- unique(all.dat$spc)
c.to.K <- function(x) x + 273.15
for (i in 1:length(spc.list)) {
  spc.index <- which(all.dat$spc == spc.list[i])
  tmp.dat <- dplyr::filter(all.dat, spc == spc.list[i])
  tmp.dat$po2 <- DO.unit.convert(tmp.dat$lc50,
                                 DO.units.in = tmp.dat$units[1],
                                 DO.units.out = "PP",
                                 bar.press = 1,
                                 bar.units.in= "atm",
                                 temp.C = tmp.dat$temp,
                                 bar.units.out = "kpa",
                                 salinity = tmp.dat$salinity[1],
                                 salinity.units = "uS")
  all.dat$po2[spc.index] <- tmp.dat$po2
}

all.dat$spc <- as.factor(all.dat$spc)
nspc <- length(levels(all.dat$spc))


tref <- 15
all.dat$inv.temp <- (1 / kb)  * ( 1 / (all.dat$temp + 273.15) - 1 / (tref + 273.15))

### Fit Models ####
fit <- glmmTMB(-log(po2) ~ diag(1+inv.temp|spc) + log(b), data = all.dat,
               family = gaussian(link = "identity"),
               se = TRUE
)

fit2 <- glmmTMB(-log(po2) ~ (1|spc) + log(b) + inv.temp, data = all.dat,
                family = gaussian(link = "identity"),
                se = TRUE
)
AICmat <- c(AIC(fit), AIC(fit2))
print(AICmat)
## model with fixed effect fits much better
# Extract parameter estimates
spc <- levels(all.dat$spc)
n <- fixef(fit2)$cond[2]
Eo <- fixef(fit2)$cond[3]
Aospc.effects <- exp(fixef(fit2)$cond[1] + ranef(fit2)$cond$spc$`(Intercept)`)

all.dat$bnpo2 <- log(all.dat$b^n * all.dat$po2)






species.names <- c("Atlantic cod", "Common eelpout", "Atlantic silverside", "N. swellfish", "Sablefish", "Sharpsnout sea bream", "Winter flounder")
for (i in 1:length(spc)) all.dat$Ao[all.dat$spc == spc[i]] <- Aospc.effects[i]
all.dat$Eo <- Eo

all.dat$bnpo2hat <- with(all.dat, -Eo * inv.temp - log(Ao) )

## TODO: Visualize fit for sablefish for range of temperatures observed in data
new.dat <- data_frame(spc = "sablefish",
                      b = 675,
                      inv.temp =(1/kb) * seq(1/c.to.K(13.6) - 1/c.to.K(15),
                                             1/c.to.K(2.98) - 1/c.to.K(15),
                                             length.out = 10)
)                                                                   

fitted <- predict(fit, newdata = new.dat)
new.dat$po2 <- exp(-fitted)
new.dat$Ao <- Aospc.effects[5]
new.dat$Eo <- Eo
new.dat$bnpo2hat <- with(new.dat, -Eo * inv.temp - log(Ao))

ggplot(data = all.dat, aes(x = inv.temp, y=(bnpo2), col = spc)) + 
  geom_point(size = 3) +
  geom_line(size = 1, aes(y = bnpo2hat)) + 
  geom_line(data = new.dat, aes(x = inv.temp, y = bnpo2hat)) +
  scale_color_brewer(palette = "Accent",name = "", labels = species.names) +
  labs(x = bquote("Inverse Temperature"~(k[b]^{-1}~ (T^-1 - T[ref]^-1))), y = bquote('log'~(B^n~pO[2]))) +
theme_bw() +
  theme(panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text= element_text(size = 10)) +
  theme(legend.position = c(0.8, 0.875))

ggsave(filename = "plots/sablefish_mi_parameters.png", dpi = 300, device = "png", height = 6, width = 6, units = "in")

