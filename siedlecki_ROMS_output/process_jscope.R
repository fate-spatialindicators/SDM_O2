# Extract J-SCOPE ROMS output from FATE spatial indicators G Drive and reformat as long data frame
# Code written by Halle Berger and Lewis Barnett

library(googledrive)

# this should fire up browser window or prompt in the console
drive_auth()

# grab data (note that file sizes are huge for anything but the annual data)
#target = drive_ls(as_id("https://drive.google.com/drive/u/0/folders/10VmLS2h2TCCVxClBqjTShm8jCuBIPxQ7")) #annual
#target = drive_ls(as_id("https://drive.google.com/drive/u/0/folders/1hXEC_-AizPaE9dGJ2SMIwJpfQdToxy5H")) #seasonal
#target = drive_ls(as_id("https://drive.google.com/drive/u/0/folders/15PJPu3GmUSQeOfAOKF12PEtcgkwUrJSC")) #monthly
target = drive_ls(as_id("https://drive.google.com/drive/u/0/folders/1IN2kPJwWXerl44vyXabg7jAE2JOFiiCM")) #daily

# pull in each file as dribble
for(i in 1:nrow(target)) {
  drive_download(file=target[i,], 
                 path=paste0("siedlecki_ROMS_output/",target$name[i]), 
                 overwrite = TRUE)
}

# read matlab files
library(R.matlab)
years <- 2009:2018 # update years as needed but only these are present in the saved outputs as of 4/5/21
#mnths <- month.abb
#dat_yr <- data.frame(lat = NA, lon = NA, bottom_temp = NA, bottom_o2 = NA, year = NA)
#dat_mnth <- data.frame(lat = NA, lon = NA, bottom_temp = NA, bottom_o2 = NA, year = NA, month = NA)
dat_day <- data.frame(lat = NA, lon = NA, bottom_temp = NA, bottom_o2 = NA, date = NA)
Matlab2Rdate <- function(val) as.Date(val - 1, origin = "0000-01-01") #need for daily files

# extract and reformat into long data frame with all time units
for(i in 1:length(years)){ #annual, monthly, & daily
  #for(j in 1:length(mnths)){ #monthly
  dat_tmp <- readMat(paste0("siedlecki_ROMS_output/",years[i],".mat")) #annual & daily
  #dat_tmp <- readMat(paste0("siedlecki_ROMS_output/",mnths[j],years[i],".mat")) #monthly
    for(j in 1:length(dat_tmp$time)){ #daily
    #coords <- expand.grid(lat = dat_tmp$coords[[1]][,1], lon = dat_tmp$coords[[2]][1,], KEEP.OUT.ATTRS = FALSE) #annual & monthly
    coords <- expand.grid(lat = dat_tmp$coords[[2]][j,,1], lon = dat_tmp$coords[[3]][j,1,], KEEP.OUT.ATTRS = FALSE) #daily
  #dat_tmp_df <- cbind(coords, bottom_temp = c(dat_tmp$bt), bottom_o2 = c(dat_tmp$bottoxy), year = years[i]) #annual
  #dat_tmp_df <- cbind(coords, bottom_temp = c(dat_tmp$bt), bottom_o2 = c(dat_tmp$bottoxy), year = years[i], month = mnths[j]) #monthly
  dat_tmp_df <- cbind(coords, bottom_temp = c(dat_tmp$bt[j,,]), bottom_o2 = c(dat_tmp$bottoxy[j,,]), date = Matlab2Rdate(dat_tmp$time[j])) #daily
  #dat_yr <- rbind(na.omit(dat_yr), na.omit(dat_tmp_df))
  #dat_mnth <- rbind(na.omit(dat_mnth), na.omit(dat_tmp_df))
  dat_day <- rbind(na.omit(dat_day), na.omit(dat_tmp_df))
  }
}

#summary(dat_yr)
#summary(dat_mnth)
summary(dat_day)
#saveRDS(dat_yr, "siedlecki_ROMS_output/annual_df.rds")
#saveRDS(dat_mnth, "siedlecki_ROMS_output/monthly_df.rds")
saveRDS(dat_day, "siedlecki_ROMS_output/daily_df.rds")
