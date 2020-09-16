load_data <- function() {
  dat = readRDS("survey_data/joined_nwfsc_data.rds")
  # analyze sablefish for years and hauls with adequate oxygen and temperature data, within range of occurrence
  dat = filter(dat, species == "sablefish", year%in%seq(2010,2015), 
               !is.na(temp), !is.na(o2), !is.na(sal),
               latitude_dd > min(latitude_dd[which(cpue_kg_km2>0)]),
               latitude_dd <= max(latitude_dd[which(cpue_kg_km2>0)]),
               longitude_dd > min(longitude_dd[which(cpue_kg_km2>0)]),
               longitude_dd < max(longitude_dd[which(cpue_kg_km2>0)]))
  
  # compute metabolic index (mi) --------------------------------------------
  # converted from Halle Berger matlab script
  
  #O2 from trawl data is in ml/l - may need to be converted to umol/kg
  gas_const = 8.31
  partial_molar_vol = 0.000032
  kelvin = 273.15
  boltz = 0.000086173324
  
  #calculate percent saturation for O2 - assumes  units of mL O2/L
  # Input:       S = Salinity (pss-78)
  #              T = Temp (deg C) ! use potential temp
  #depth is in meters
  #[umole/kg] = [ml/L]*44660/(sigmatheta(P=0,theta,S) + 1000)
  dat$SA = gsw_SA_from_SP(dat$sal,dat$depth,dat$longitude_dd,dat$latitude_dd) #absolute salinity for pot T calc
  dat$pt = gsw_pt_from_t(dat$SA,dat$temp,dat$depth) #potential temp at a particular depth
  dat$CT = gsw_CT_from_t(dat$SA,dat$temp,dat$depth) #conservative temp
  dat$sigma0 = gsw_sigma0(dat$SA,dat$CT)
  dat$o2_umolkg = dat$o2*44660/(dat$sigma0+1000)
  # calc o2 solubility, relies on o2 in umol/kg
  gsw_O2sol_SP_pt <- function(sal,pt) {
    x = dat$sal
    pt68 = dat$pt*1.00024
    y = log((298.15 - pt68)/(273.15 + pt68))
    
    a0 =  5.80871
    a1 =  3.20291
    a2 =  4.17887
    a3 =  5.10006
    a4 = -9.86643e-2
    a5 =  3.80369
    b0 = -7.01577e-3
    b1 = -7.70028e-3
    b2 = -1.13864e-2
    b3 = -9.51519e-3
    c0 = -2.75915e-7
    
    O2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
    return(O2sol)
  }
  
  dat$O2_Sat0 = gsw_O2sol_SP_pt(dat$sal,dat$pt)
  
  #= o2satv2a(sal,pt) #uses practical salinity and potential temp - solubity at p =1 atm
  dat$press = exp(dat$depth*10000*partial_molar_vol/gas_const/(dat$temp+kelvin))
  dat$O2_satdepth = dat$O2_Sat0*dat$press
  
  #solubility at p=0
  dat$sol0 = dat$O2_Sat0/0.209
  dat$sol_Dep = dat$sol0*dat$press
  dat$po2 = dat$o2_umolkg/dat$sol_Dep
  
  # species-specific parameters
  Ao = 1.16625e-13
  Eo = 0.8736 * 0.5 # from cod, 0.8736.  Make it one half or double
  B = 3000 # size in grams, roughly average (initial calculations used 10000g)
  n = -0.208 # borrowed from cod 
  
  dat$mi = B^n*Ao*dat$po2/exp(-1*Eo/(boltz*(dat$temp+kelvin)))
  
  # prepare data and models -------------------------------------------------
  
  dat <- select(dat, species, year, longitude_dd, latitude_dd, cpue_kg_km2,
                o2, temp, depth, mi, po2)
  
  
  # UTM transformation
  dat_ll = dat
  coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
  proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = spTransform(dat_ll, 
                        CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat = as.data.frame(dat_utm)
  dat = dplyr::rename(dat, longitude = longitude_dd, 
                      latitude = latitude_dd)
  return(dat)
}

back.convert <- function(x, mean_orig, sd_orig) (x* sd_orig+mean_orig)

plot_map <- function(dat, column = "est") {
  ggplot(dat, aes_string("longitude", "latitude", fill = column)) +
    geom_tile() +
    facet_wrap(~year) +
    coord_fixed()
}

get_models <- function() {
  formula <- c("depth + I(depth^2) + as.factor(year)", 
                                      "depth + I(depth^2) + as.factor(year) + temp", 
                                      "depth + I(depth^2) + as.factor(year) + o2",
                                      "depth + I(depth^2) + as.factor(year) + po2",
                                      "depth + I(depth^2) + as.factor(year) + mi",
                                      "depth + I(depth^2) + as.factor(year) + temp + o2",
                                      "depth + I(depth^2) + as.factor(year) + temp + o2 + temp*o2",
                                      "depth + I(depth^2) + as.factor(year) + temp + po2",
                                      "depth + I(depth^2) + as.factor(year) + temp + po2 + temp * po2",
                                      "depth + I(depth^2) + as.factor(year) + breakpt(o2)",
                                      "depth + I(depth^2) + as.factor(year) + breakpt(o2)+ temp",
                                      "depth + I(depth^2) + as.factor(year) + breakpt(po2)",
                                      "depth + I(depth^2) + as.factor(year) + breakpt(po2) + temp",
                                      "depth + I(depth^2) + as.factor(year) + breakpt(mi)"
                                      )
}
                      