#Xtracto functions to extract ROMS data
#Functions initially built by Jordan Watson with Jameal Samhouri cheerleading

library(tidync)
library(tidyverse)
library(here)

# bottom O2
nc <- tidync(here::here(
  "jacox_ROMS_output",
  "oxygen_bottom_era5_monthly_1980_2010.nc"
))

#  If we look at nc, we see that the "active grid" is the one that gives us only index values.
#  This seems like a space saving approach - more compact to use index values.
nc

#  So let's see if the grid that only contains lon,lat,depth includes the values we want.
#  To do this we have to activate one of the other grids.
mygrid <- nc %>%
  activate("D0,D1") %>%
  hyper_tibble()

#  Cool - looks like our lat-longs were hiding in a different grid.
head(mygrid)

#  Ha! I bet they did the same thing with time. Let's look at the time grid.
#  Faster to create the POSIXct field here
mytime <- nc %>%
  activate("D2") %>%
  hyper_tibble() %>%
  mutate(xdate=as.POSIXct(paste(year,month,"01",sep="-"),format="%Y-%m-%d",tz="UTC"))

head(mytime)

#  Now let's create our tibble that includes our variable of interest.
data <- nc %>%
  hyper_tibble(select_var="oxygen")

#  We should just be able to join these all together.
test <- data %>%
  inner_join(mygrid) %>%
  inner_join(mytime)

#  Or in one piped action...
nc %>%
  hyper_tibble(select_var="oxygen") %>%
  inner_join(nc %>%
               activate("D0,D1") %>%
               hyper_tibble()) %>%
  inner_join(nc %>%
               activate("D2") %>%
               hyper_tibble())

#  Or in one piped action...
#  Make sure that the as.POSIXct is inside the inner_join instead of outside or it'll take forever.
nc %>%
  hyper_tibble(select_var="oxygen") %>%
  inner_join(nc %>%
               activate("D0,D1") %>%
               hyper_tibble()) %>%
  inner_join(nc %>%
               activate("D2") %>%
               hyper_tibble() %>%
               mutate(xdate=as.POSIXct(paste(year,month,"01",sep="-"),format="%Y-%m-%d",tz="UTC")))