# get ROMS bottom temp, convert to flat df
# get joined_nwfsc_data 
# make sure coordinates and dates make sense for both df's
  # add month_year to joined_nwfsc_data df
# join ROMS bottom temp data to joined_nwfsc_data by month_year

library(ncdf4)

# get ROMS bottom temp, convert to flat df

wc_temp_bottom_file <-  here::here(
  "jacox_ROMS_output",
  "temp_bottom_era5_monthly_1980_2010.nc"
)


ThisNC.nc<-nc_open(wc_temp_bottom_file) 

save(ThisNC.nc,
     file=here::here(
       "jacox_ROMS_output",
       "temp_bottom_era5_monthly_1980_2010.Rdata")
)

# could start here on 11 feb 2020. read in .rdata instead of netcdf. then turn it into a flat file with dd lat long year month temp. then join

# JUST LIST ALL THE VARIABLES IN THE NCDF FILE
# THIS CODE JUST PRINTS VARIABLE NAMES TO SCREEN, BUT CAN BE USEFUL FOR MAKING THE .CSV INPUT FILES. 
print('LISTING ALL VARIABLE NAMES IN THE NETCDF FILE')
numVarsInNC<-length(ThisNC.nc$var)
numVarsInNC
ListOfVarsFromNC<-matrix(numVarsInNC,1)
ListOfVarsFromNC

for (groupIndex in 1:numVarsInNC )   
{
  thisVar <- ThisNC.nc$var[[groupIndex]]  
  #print(thisVar$name)    
  ListOfVarsFromNC[groupIndex]<-thisVar$name
}
write.csv(ListOfVarsFromNC,file=here::here(
  "jacox_ROMS_output",
  "ListOfVarsFromNC.csv")
)
print('THE LISTING ABOVE IS ALL VARIABLE NAMES IN CDF')
print('THIS MAY BE USEFUL FOR MAKING .CSV INPUT FILES')
print(' THE LIST OF VARIABLES IS SAVED AS ListOfVarsFromNC.csv')

# variables are: 
# "temp"  "lon"   "lat"   "depth" "year"  "month"

lat <- ncvar_get( ThisNC.nc,"lat") # extract the data from the variable. The variable contains lots of other metainfo like units, name, etc.
latDims<-dim(lat)  # Just use volume to see how many time steps are in the data

lon<- ncvar_get( ThisNC.nc,"lon") # extract the data from the variable. The variable contains lots of other metainfo like units, name, etc.
lonDims<-dim(lon)  # Just use volume to see how many time steps are in the data

temp_bottom <- ncvar_get( ThisNC.nc,"temp") # extract the data from the variable. The variable contains lots of other metainfo like units, name, etc.
temp_bottomDims<-dim(temp_bottom)  # Just use volume to see how many time steps are in the data

year<- ncvar_get( ThisNC.nc,"year") # extract the data from the variable. The variable contains lots of other metainfo like units, name, etc.
yearDims<-dim(year)  # Just use volume to see how many time steps are in the data


month<- ncvar_get( ThisNC.nc, "month" ) # extract the data from the variable. The variable contains lots of other metainfo like units, name, etc.
monthDims<-dim(month)  # Just use volume to see how many time steps are in the data

source(here::here(
  "jacox_ROMS_output",
  "PlotNCmap.r")
)
PlotNCmap(temp_bottom[,,6],"temp_bottom",lon,lat)






#library(sm)  # for pause()


#library(chron) #added 4/30/12
#remove.packages("compositions") # for mean.row
#remove.packages("fields") # for image.plot()
#library(graphics)
#remove.packages("MASS") # for write.matrix()
#remove.packages("matlab")
#library(grDevices)

#---------------------
# USER CAN DEFINE INPUT FILES, STARTYEAR, TIMESTEP, SAVE DIRECTORY, ETC. HERE:  
#Here are some things that are hardcoded, but change here if you want. 
#boxesToUse<-7

nc_close(ThisNC.nc)



##############################################
##############################################
###### GRAVEYARD
##############################################
##############################################

## jameal tried to implement tidync and failed. for now.
library(tidync)

# example
# file <- system.file("extdata", "oceandata", "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc", package = "tidync")
# library(tidync)
# tidync(file) 
# browseVignettes(package = "tidync")

wc_temp_bottom_file <-  here::here(
  "jacox_ROMS_output",
    "temp_bottom_era5_monthly_1980_2010.nc"
  )
tidync(wc_temp_bottom_file)
tidync(wc_temp_bottom_file) %>% hyper_dims()
tidync(wc_temp_bottom_file) %>% hyper_grids()
tidync(wc_temp_bottom_file) %>% activate(temp)

str(tidync(wc_temp_bottom_file) %>% hyper_filter())

lat <- ncvar_get( tidync(wc_temp_bottom_file) ,"lat")

tidync(wc_temp_bottom_file) %>%
  hyper_filter(latitude = between(latitude, -78, -75.8), 
               longitude = between(longitude, 165, 171))


(subs <- tidync(wc_temp_bottom_file) %>% hyper_filter(longitude = longitude < 126, latitude = latitude > 46, time = time > 300))
subs %>% hyper_tibble()

##############################################
##############################################
##############################################
