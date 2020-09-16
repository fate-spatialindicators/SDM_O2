library(googledrive)

# this should fire up browser window or prompt in the console
drive_auth()

# here grabbing annual data to play with smaller file size, other file sizes are huge
target = drive_ls(as_id("https://drive.google.com/drive/u/0/folders/10VmLS2h2TCCVxClBqjTShm8jCuBIPxQ7")) #annual
#target = drive_ls(as_id("https://drive.google.com/drive/u/0/folders/1hXEC_-AizPaE9dGJ2SMIwJpfQdToxy5H")) #seasonal
#target = drive_ls(as_id("https://drive.google.com/drive/u/0/folders/15PJPu3GmUSQeOfAOKF12PEtcgkwUrJSC")) #monthly
#target = drive_ls(as_id("https://drive.google.com/drive/u/0/folders/1IN2kPJwWXerl44vyXabg7jAE2JOFiiCM")) #daily

# pull in each file as dribble
for(i in 1:nrow(target)) {
  drive_download(file=target[i,], 
                 path=paste0("siedlecki_ROMS_output/",target$name[i]), 
                 overwrite = TRUE)
}

# read matlab files
library(R.matlab)
years <- 2009:2018
dat_yr <- data.frame(lat = NA, lon = NA, bottom_temp = NA, bottom_o2 = NA, year = NA)

# extract and reformat into long data frame with all time units
for(i in 1:length(years)){
  dat_tmp <- readMat(paste0("siedlecki_ROMS_output/",years[i],".mat"))
  #str(dat_tmp) # output is list of matrices associated with $coords
  #str(dat_tmp$coords[[2]])
  dat_tmp$coords[[1]][,1]#unique lat
  dat_tmp$coords[[2]][1,]#unique lon
  coords <- expand.grid(lat = dat_tmp$coords[[1]][,1], lon = dat_tmp$coords[[2]][1,], 
                      KEEP.OUT.ATTRS = FALSE)
  dat_tmp_df <- cbind(coords, bottom_temp = c(dat_tmp$bt), bottom_o2 = c(dat_tmp$bottoxy), year = years[i])
  dat_yr <- rbind(na.omit(dat_yr), na.omit(dat_tmp_df))
}

summary(dat_yr)
saveRDS(dat_yr, "siedlecki_ROMS_output/annual_df.rds")
