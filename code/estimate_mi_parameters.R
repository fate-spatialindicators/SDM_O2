library(ggplot2)
library(viridis)
library(gsw)
library(dplyr)
library(rMR)
library(lme4)
library(RColorBrewer)
library(rstan)
library(shinystan)

library(rMR)
kb <-  8.617333262145E-5
all.dat<- read.csv("data/metabolic_index_data.csv", header = T)
all.dat$po2 <- NA
# cycle by species and measurement type, and get po2 for each
spc.list <- unique(all.dat$spc)

for (i in 1:length(spc.list)) {
  spc.index <- which(all.dat$spc == spc.list[i])
  tmp.dat <- dplyr::filter(all.dat, spc == spc.list[i])
  tmp.dat$po2 <- DO.unit.convert(tmp.dat$lc50,
                                 DO.units.in = tmp.dat$units[1],
                                 DO.units.out = "PP",
                                 bar.press = 1,
                                 bar.units.in= "atm",
                                 temp.C = tmp.dat$temp,
                                 bar.units.out = "atm",
                                 salinity = tmp.dat$salinity[1],
                                 salinity.units = "uS")
  all.dat$po2[spc.index] <- tmp.dat$po2
}

all.dat$spc <- as.factor(all.dat$spc)
# fit global model 
tref <- 15
all.dat$inv.temp <- (1 / kb)  * ( 1 / (all.dat$temp + 273.15) + 1 / (15 + 273.15))


all.dat$spc <- as.factor(all.dat$spc)
fit <- lmer(log(po2) ~inv.temp + log(b) + (1|spc), data = all.dat)


Eo <- - fixef(fit)[2]
n <- - fixef(fit)[3]

spc.effects <- coef(fit)


species.names <- c("Atlantic cod", "Common eelpout", "Atlantic silverside", "N. swellfish", "Sablefish", "Sharpsnout sea bream", "Winter flounder")
all.dat$bnpo2 <- all.dat$b^n * all.dat$po2

ggplot(data = all.dat, aes(x = inv.temp, y=log(bnpo2), col = spc)) + 
  geom_point(size = 3) +
  scale_color_brewer(palette = "Accent",name = "", labels = species.names) +
  labs(x = bquote("Inverse Temperature"~(k[b]~T)^-1), y = bquote('log'~(B^n~pO[2]))) +
theme_bw() +
  theme(panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) +
  theme(legend.text= element_text(size = 12)) +
  theme(legend.position = c(0.8, 0.875))

ggsave(filename = "plots/sablefish_mi_parameteris.png", dpi = 300, device = "png", height = 6, width = 6, units = "in")


## Repeat using bayesian Model
spc <- as.numeric(all.dat$spc)
inv_temp <- all.dat$inv.temp
logb <- log(all.dat$b)
y <- log(all.dat$po2)
N <- length(y)
n <- length(unique(spc))

stan.model <- "Code/Stan/fit_mi.stan"
stan.pars <- c("Eo", "n", "Ao", "mu_Ao", 'logsigma_Ao', 'sigma_Ao', 'sigma')
stan.data <- list(N = N, nspc = n, y = y, logb = logb, inv_temp = inv_temp, spc = spc)
n_chains <- 3

init.chain <- function(chain_id = 1, n) {
  list(
    Eo = rnorm(1, 0, 1),
    b = rnorm(1,0,1),
    Ao_raw = rnorm(n, 0, 1),
    mu_Ao = rnorm(1, 30, 20),
    logsigma_Ao = rnorm(1, 0, 1),
    sigma = runif(1, 0.1, 10)
  )
}


init_ll <- lapply(1:n_chains, function(id)
  init.chain(chain_id = id, n = n))



rstan_options(auto_write = TRUE)  # this option stops Stan from re-compiling if not necessary
options(mc.cores = parallel::detectCores()) # this is nice because it allows each chain to be run in parallel on a separate core of your processor

niters <- 5000 # how long should each chain be.  1000 is probably fine

fit <- stan(
  file = stan.model,
  data = stan.data,
  iter = niters,
  pars = stan.pars,
  warmup = floor(niters / 2),
  chains = n_chains,
  thin = 5,
  algorithm = 'NUTS',
  init = init_ll,
  verbose = FALSE,
  control = list(adapt_engaged = TRUE, adapt_delta = 0.9, max_treedepth = 15)
)

params <- extract(fit)

main.params <- cbind(params$Eo, params$n, params$mu_Ao, params$logsigma_Ao)

param_means <- colMeans(main.params)
param_sigma <- cov(main.params)
save(file = "mi_parameters_bayesian.Rdata", fit)

#launch_shinystan(fit)
