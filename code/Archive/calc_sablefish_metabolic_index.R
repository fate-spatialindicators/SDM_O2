library(rMR)

B.sable <- 0.675
o2sat.crit <- 5.4 # percent saturation at loss of equilibrium
temp <- 12 # Paper claims this is the holding temperature (between 10 - 11 degrees).
sal <-28 # rough guess, if water coming from Puget Sound, it'll be around here.

po2.crit.sable <- DO.unit.convert(o2sat.crit, 
                            DO.units.in = "pct",
                            DO.units.out = "PP",
                            bar.press = 1,
                            bar.units.in= "atm",
                            temp.C = temp,
                            bar.units.out = "atm",
                            salinity = sal,
                            salinity.units = "uS")
po2.threshold <- DO.unit.convert(24, 
                                 DO.units.in = "pct",
                                 DO.units.out = "PP",
                                 bar.units.in = "atm",
                                 bar.press = 1,
                                 temp.C = temp,
                                 bar.units.out = "atm",
                                 salinity = sal)
## Get ready to plot
## function to plot phi for Cod
kb <- 8.617333262145E-5
Eo <- 0.4903464 
Ao <- 1.493868e-07 

B <- 1000
n <- -0.3512716

t.range <- 5:17
t.ref <- 15
k.range <- t.range + 273.15
k.ref <- t.ref + 273.15

inv.k <- 1/k.range
log.po2.crit <- -log(Ao) -Eo/kb * (1/k.range) # Deutsch et al


inverse.temp <- 1/ (kb * k.range)

plot(inverse.temp, log.po2.crit,
     type = "l",
     lwd = 2,
     ylab = "log pO2 crit",
     xlab = "inverse T",
     ylim = c(-6,-2),
     xlim = c(39.5, 42.5),
     las =1,
     col = "red")

cod.data <- matrix(c(40.000,	-3.785,
                     40.270,	-4.365,
                     40.837,	-4.413,
                     40.981,	-4.522,
                     41.128,	-4.762,
                     41.714,	-5.542), 
                   byrow = T,
                   nrow = 6, ncol = 2)
points(cod.data[,1], cod.data[,2],
       pch = 21,
       bg = "red")

points(1/((temp + 273.15) * kb), log(po2.crit.sable * B.sable^n),
          pch = 21,
          bg = "purple")
# calculate Ao that makes line go exactly through the single data point
newEo <- 1 * Eo
newlogAo <- -log(po2.crit.sable * B.sable^n) -newEo* 1/((temp + 273.15) * kb)
newAo <- exp(newlogAo)
lines(inverse.temp, -newlogAo-newEo*inverse.temp, lwd = 2, col = "purple")

