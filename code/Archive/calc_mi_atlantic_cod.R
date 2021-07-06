# re-create atlantic cod estimates
kb <-  8.617333262145E-5
plante.df <- data.frame(b = c(1821,1653,622, 517,1734,1842,697,880),
                        temp = c(rep(6,4), rep(2,4)),
                        lc50 = c(26.4, 22.3, 18.7, 19.2, 22, 20, 23.7, 21)
)
ss.df <- data.frame(b = rep(145, times = 6),
                    temp = c(17, 15, 11, 10, 9, 5),
                    lc50 = c(27, 14, 13, 12, 10, 6)
)

cod.df <- rbind(plante.df, ss.df)

library(rMR)
cod.df$po2 <- DO.unit.convert(cod.df$lc50,
                                 DO.units.in = "pct",
                                 DO.units.out = "PP",
                                 bar.press = 1,
                                 bar.units.in= "atm",
                                 temp.C = cod.df$temp,
                                 bar.units.out = "atm",
                                 salinity = 33,
                                 salinity.units = "uS")

cod.df$inv.temp <- 1 / (kb * (cod.df$temp+273.15))

fit <- lm(log(po2) ~  inv.temp + log(b), data = cod.df)
coef(fit)
Ao <-exp(-12.4282)
Eo <- 0.449
n <- -0.465

# MI for Sablefish data
mi <- .01 * Ao * .675^n * exp(Eo / (kb * (12 + 273.5)))
Ao.sable <-  Ao / mi.sable
mi.sable <- .01 * Ao.sable * .675^n * exp(Eo / (kb * (12 + 273.5)))
cod.df$o2 <- DO.unit.convert(cod.df$lc50,
                             DO.units.in = "pct",
                             DO.units.out = "mg/L",
                             bar.press = 1,
                             bar.units.in= "atm",
                             temp.C = cod.df$temp,
                             salinity = 33,
                             salinity.units = "uS")

ss.df$po2 <- DO.unit.convert(ss.df$lc50,
                             DO.units.in = "pct",
                             DO.units.out = "PP",
                             bar.press = 1,
                             bar.units.in= "atm",
                             temp.C = ss.df$temp,
                             bar.units.out = "atm",
                             salinity = 35,
                             salinity.units = "uS")
ss.df$o2 <- DO.unit.convert(ss.df$lc50,
                             DO.units.in = "pct",
                             DO.units.out = "mg/L",
                             bar.press = 1,
                             bar.units.in= "atm",
                             temp.C = ss.df$temp,
                             salinity = 35,
                             salinity.units = "uS")
ss.df$inv.temp <- 1 / (kb * (ss.df$temp+273.15))
fit <- lm(log(po2) ~inv.temp, data = ss.df)
coef(fit)
