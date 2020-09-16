#Xtracto functions to extract ROMS data
#Functions initially built by K.Scales and modified by S.Brodie,  H.Welch, and E.Hazen

#Xtracto_ROMs mean and sd function
#Set desired resolution in which to apply a function (mean or sd), with 'name' indicating the text string of the function
# EXAMPLE: getvarROMS(ROMS_files_hist[1], 'curl', input_file, desired.resolution = 0.3, mean, 'mean_0.3')

inpts=read.csv("/Users/heatherwelch/Downloads/jam.csv")
nc="/Users/heatherwelch/Downloads/temp_bottom_era5_monthly_1980_2010.nc"
varname="temp"
name="temperature"

getvarROMS <- function(nc,varname,inpts,name){
  inpts=inpts %>% mutate(dt=substr(date,1,8)) %>% mutate(dt=paste0(dt,"01"))
  inpts$dt <- as.POSIXct(inpts$dt, '%Y-%m-%d',tz='UTC')
  nc.data <- nc_open(nc, write=FALSE)
  lat <- ncvar_get(nc.data,'lat'); lat <- lat[1,]
  lon <- ncvar_get(nc.data,'lon'); lon <- lon[,1]
  nrows <- length(lat); ncols <- length(lon)
  # mth <- ncvar_get(nc.data,'month')
  yr <- ncvar_get(nc.data,'year'); mth <- ncvar_get(nc.data,'month');day=rep("01",length(mth))
  tim <- as.POSIXct(paste(yr,mth,day,sep='-'),tz='UTC')
  for (i in 1:nrow(inpts)){
    print(paste(varname,inpts$date[i],sep=' '))
    if (inpts$dt[i] %in% tim){
      xdate <- which(inpts$dt[i]==tim)
      c <- which.min(abs(lon-inpts$lon[i]))
      r <- which.min(abs(lat-inpts$lat[i]))
        data.var.point  <-  ncvar_get(nc.data,varname,start=c(c,r,xdate),
                                      count=c(1,1,1),verbose=FALSE)
        inpts[i,paste(varname,'_',name,sep='')] <- data.var.point
     
    }
  }
  nc_close(nc.data)
  return(inpts)
}

test=getvarROMS(nc=nc,varname = varname,inpts = inpts,name=name)
