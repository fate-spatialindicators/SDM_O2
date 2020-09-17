devtools::install_github("pbs-assess/sdmTMB")
library(sdmTMB)
library(dplyr)
library(sp)
library(gsw)

env_data <- readRDS("survey_data/bc-synoptic-env.rds")
trawl_data <- readRDS("survey_data/bc-synoptic-trawls.rds")
dat_all <- left_join(trawl_data, env_data, by = "fishing_event_id")

# renaming to match wc data
dat_all <- rename(dat_all, o2 = do, temp = temperature, cpue_kg_km2 = density_kgpm2, sal = salinity, 
              depth = depth_m, longitude_dd = longitude, latitude_dd = latitude)

sp_survey <- arrange(expand.grid(sp = c("sablefish","pacific cod"), survey = unique(dat_all$survey)), sp)

for(j in 1:nrow(sp_survey)){
  print(paste0("model set # ", j, " of ", nrow(sp_survey)))
  # analyze years and hauls with adequate oxygen and temperature data, within range of occurrence
  dat = filter(dat_all, species == sp_survey$sp[j], survey == sp_survey$survey[j],
               !is.na(temp), !is.na(o2), !is.na(sal), !is.na(depth),
               latitude_dd > min(latitude_dd[which(cpue_kg_km2>0)]),
               latitude_dd <= max(latitude_dd[which(cpue_kg_km2>0)]),
               longitude_dd > min(longitude_dd[which(cpue_kg_km2>0)]),
               longitude_dd < max(longitude_dd[which(cpue_kg_km2>0)]))
  
  #check on sampling coverage by survey
  #print(dat %>% count(survey, year), n = 25)
  
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
  if(sp_survey$sp[j] == "sablefish"){
    Ao = 1.16625e-13 # sablefish specific
    B = 1200 # size in grams, average in survey
    Eo = 0.8736 * 0.5 # from cod, 0.8736.  Make it one half or double
  }
  else{
    Ao = 3.11E-14 # cod
    B = 1250 # eyeball estimate from Anderson et al. 2019 report (varies by survey in reality)
    Eo = 0.8736
  }
  n = -0.208 # borrowed from cod 
  
  dat$mi = B^n*Ao*dat$po2/exp(-1*Eo/(boltz*(dat$temp+kelvin)))
  
  # prepare data and models -------------------------------------------------
  
  dat <- select(dat, species, year, longitude_dd, latitude_dd, cpue_kg_km2,
                o2, temp, depth, mi, po2)
  
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
  
  # UTM transformation
  dat_ll = dat
  coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
  proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = spTransform(dat_ll, 
                        CRS("+proj=utm +zone=9 +datum=WGS84 +units=km")) # check on zone
  # convert back from sp object to data frame
  dat = as.data.frame(dat_utm)
  dat = dplyr::rename(dat, longitude = longitude_dd, 
                      latitude = latitude_dd)
  
  # create combination of covariates and threshold responses for different models
  m_df = data.frame(
    spatial_only = rep(TRUE,14), 
    depth_effect = rep(TRUE,14),
    time_varying = rep(FALSE,14),
    threshold_function = c(rep("NA",9),rep("linear", 5)),
    covariate1 = c("none","temp","o2","po2","mi",rep("temp",4),rep("o2",2), rep("po2", 2),"mi"),
    covariate2 = c(rep("none",5),"o2","o2","po2", "po2","none", "temp", "none", "temp","none"),
    interaction = c(rep(FALSE,6),TRUE, FALSE, TRUE, rep(FALSE,5)),
    tweedie_dens = rep(NA,14) # set up vector to store performance data if using cv
  )
  
  # run models for each combination of settings/covariates in df ------------
  
  use_cv = FALSE # specify whether to do cross validation or not
  spde <- make_spde(x = dat$longitude, y = dat$latitude, n_knots = 250) # choose # knots
  
  # fit models with typical approach
  for(i in 1:nrow(m_df)){
    print(paste0("model # ", i, " of ", nrow(m_df)))
    
    formula = paste0("cpue_kg_km2 ~ 0 + as.factor(year)")
    
    # format data and formula based on combination of arguments in model settings df
    if(m_df$threshold_function[i] == "linear") {
      formula = paste0(formula, " + ", "breakpt(enviro1)")
    }
    if(m_df$threshold_function[i] == "logistic") {
      formula = paste0(formula, " + ", "logistic(enviro1)")
    } 
    if(m_df$threshold_function[i] == "quadratic") {
      formula = paste0(formula, " + ", "I(enviro1^2)")
    }
    if(m_df$threshold_function[i] == "NA") {
      formula = paste0(formula, " + ", "enviro1")
    }
      
    # add 2nd covariate if included
    if(m_df$covariate2[i] != "none") {
      # rename variables to make code generic
      sub <- dplyr::rename(dat, enviro2 = as.character(m_df$covariate2[i]))
      if(m_df$interaction[i] == TRUE){
        formula = paste0(formula, " + ", "enviro2", " + ", "enviro1", " * ", "enviro2")
      } else {
        formula = paste0(formula, " + ", "enviro2")
      }
    }
    
    if(m_df$covariate1[i] == "none") {
      formula = paste0("cpue_kg_km2 ~ 0 + as.factor(year)")
      sub <- dat
    }else{
      sub <- dplyr::rename(dat, enviro1 = as.character(m_df$covariate1[i]))
    }
    
    if(m_df$depth_effect[i]==TRUE) {
      formula = paste0(formula, " + depth + I(depth^2)")
    }
    
    # fit model with or without cross-validation
    if(use_cv==TRUE) {
      m <- try(sdmTMB_cv(
        formula = as.formula(formula),
        data = sub,
        x = "longitude", 
        y = "latitude",
        time = NULL,
        k_folds = 4,
        n_knots = 250,
        seed = 10,
        family = tweedie(link = "log"),
        anisotropy = TRUE,
        spatial_only = m_df$spatial_only[i]
      ), silent = TRUE)
      
      if(class(m)!="try-error") {
        dir_name = paste0("output/bc/",sub("pacific ","",sp_survey[j,1]),"_",sub("SYN ","",sp_survey[j,2]))
        if(!dir.exists(dir_name)){dir.create(dir_name)}
        saveRDS(m, file = paste0(dir_name,"/model_",i,"_MI_cv.rds"))
        m_df$tweedie_dens[i] = m$sum_loglik
      }
      
    } else {
      m <- try(sdmTMB(
        formula = as.formula(formula),
        data = sub,
        time = NULL,
        reml = TRUE,
        spde = spde,
        family = tweedie(link = "log"),
        anisotropy = TRUE,
        spatial_only = m_df$spatial_only[i]
      ), silent = TRUE)
      
      if(class(m)!="try-error") {
        dir_name = paste0("output/bc/",sub("pacific ","",sp_survey[j,1]),"_",sub("SYN ","",sp_survey[j,2]))
        if(!dir.exists(dir_name)){dir.create(dir_name)}
        saveRDS(m, file = paste0(dir_name,"/model_",i,"_MI.rds"))
      }
      
    }
  }
  
  # Get AIC table
  AICmat <- matrix(NA, nrow = nrow(m_df), ncol = 2)
  for (i in 1:nrow(m_df)) {
    if(file.exists(paste0(dir_name,"/model_",i,"_MI.rds"))){
      m <- readRDS(paste0(dir_name,"/model_",i,"_MI.rds"))
      AICmat[i,1] <- AIC(m)
    }else{
      AICmat[i,1] <- NA
    }
  }
  
  dAIC <- AICmat[,1] - min(AICmat[,1], na.rm = TRUE)
  dAIC <- matrix(as.numeric(sprintf(dAIC,fmt = '%.2f')), nrow = nrow(m_df), ncol = 1)
  rownames(dAIC) <- c("space + depth + year",
                      "space + depth + year + temp", 
                      "space + depth + year + o2",
                      "space + depth + year + p02",
                      "space + depth + year + mi",
                      "space + depth + year + temp + o2",
                      "space + depth + year + temp:o2",
                      "space + depth + year + temp + po2",
                      "space + depth + year + temp:po2",
                      "space + depth + year + o2(breakpoint)",
                      "space + depth + year + o2(breakpoint) + temp",
                      "space + depth + year + po2(breakpoint)",
                      "space + depth + year + po2(breakpoint) + temp",
                      "space + depth + year + mi(breakpoint)")
  sink("output/bc/aic_tables.txt", append = TRUE)
  print(sp_survey[j,])
  print(dAIC)
  sink()
}