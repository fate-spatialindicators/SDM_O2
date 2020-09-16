#Xtracto functions to extract ROMS data
#Functions initially built by K.Scales and modified by S.Brodie,  H.Welch, and E.Hazen

library(tidync)

#Xtracto_ROMs mean and sd function
#Set desired resolution in which to apply a function (mean or sd), with 'name' indicating the text string of the function
# EXAMPLE: getvarROMS(ROMS_files_hist[1], 'curl', input_file, desired.resolution = 0.3, mean, 'mean_0.3')

#inpts=read.csv("10 rows of wc survey data.csv")
# get survey data
inpts <- readRDS("survey_data/joined_nwfsc_data.rds")
inpts$date <- ymd(dat$date)
inpts$month <- month(dat$date)

# bottom temp
# nc=here::here(
#   "jacox_ROMS_output",
#   "temp_bottom_era5_monthly_1980_2010.nc"
# )
# varname="temp"
# name="ROMS_temp_bottom_era5_monthly"

# bottom O2
nc=here::here(
  "jacox_ROMS_output",
  "oxygen_bottom_era5_monthly_1980_2010.nc"
)
tidync(nc) # use this to determine varname

#tidync(nc) %>% hyper_tibble(select_var = "oxygen")
#tidync(nc) %>% hyper_tibble(select_var = "lon")

################
################
varname="oxygen" #### need to change for O2
################
################
name="ROMS_oxygen_bottom_era5_monthly"

getvarROMS <- function(nc,varname,inpts,name){
  inpts=inpts %>% mutate(dt=substr(date,1,8)) %>% mutate(dt=paste0(dt,"01")) # this turns every date into the first of the month. might want to change if i want to days between, eg june 17 and july 14 to be july but dates between july 15 and aug 14 to be aug
  inpts$dt <- as.POSIXct(inpts$dt, '%Y-%m-%d',tz='UTC')
  nc.data <- nc_open(nc, write=FALSE)
  lat <- ncvar_get(nc.data,'lat'); lat <- lat[1,]
  lon <- ncvar_get(nc.data,'lon'); lon <- lon[,1]
  nrows <- length(lat); ncols <- length(lon)
  # mth <- ncvar_get(nc.data,'month')
  yr <- ncvar_get(nc.data,'year'); mth <- ncvar_get(nc.data,'month');day=rep("01",length(mth))
  tim <- as.POSIXct(paste(yr,mth,day,sep='-'),tz='UTC')
  for (i in 1:nrow(inpts)){
    print(paste(varname,inpts$date[i],sep=' ')) # just tracking the date of each haul to match to each var
    if (inpts$dt[i] %in% tim){
      xdate <- which(inpts$dt[i]==tim) # find the ROMS data corresponding to the month of the haul
      c <- which.min(abs(lon-inpts$lon[i])) # find the closest longitude to the haul
      r <- which.min(abs(lat-inpts$lat[i])) # find the closest latitude to the haul
        data.var.point  <-  ncvar_get(nc.data,varname,start=c(c,r,xdate),
                                      count=c(1,1,1),verbose=FALSE)
        inpts[i,name] <- data.var.point # paste(varname,'_',name,sep='')
     
    }
  }
  nc_close(nc.data)
  return(inpts)
}

#test=getvarROMS(nc=nc,varname = varname,inpts = inpts,name=name)

start.time <- Sys.time()

out <- getvarROMS(nc=nc,varname = varname,inpts = inpts,name=name)

Sys.time() - start.time
# started temp_bottom ~2153 on 02112020, finished before 2236
# bottom oxygen: Time difference of 37.74415 mins

saveRDS(out, paste0("survey_data/joined_nwfsc_data",name,".rds"))
