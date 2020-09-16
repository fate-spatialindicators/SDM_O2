MatchJSCOPEnetcdfToLogbook <- function(oceanMatrixAllMonthsYears,WDFW_Crab_Logbooks_Data_ThisYearOnly,startYearForFishingSeason)
{
#------
    # GOAL IS TO match the ocean points up to logbook locations, and assign ocean variables as new columsn of the logbook rows.
   # Arguments:
  #       oceanMatrixAllMonthsYears: very large matrix of data for every grid point of ROMS Jscope NC files
  #       WDFW_Crab_Logbooks_Data_ThisYearOnly: single year of washington logbook data
  #       startYearForFishingSeason:  would be for instance 2017 for fishing year autum 2107--spring 2018.
  #                                   Important because this is used to find NC file we want to match to this logbook record.
  #
  #   Returns: WDFWThisYearLbsAndCountsConvertedToLbs,  each row has logbook record (including cath rate per trap) and associated ocean conditiosn predicted by ROMS nearest point for that month.

#----------
print("STARTING TO MATCH JSCOPE netcdf FILE TO CRAB LOGBOOK FOR 1 YEAR")

source("PlotCoastForSardine.r")

rowsForThisYear <- which(oceanMatrixAllMonthsYears[,"yearROMSstart"]==startYearForFishingSeason)


oceanMatrixThisYear <- oceanMatrixAllMonthsYears[rowsForThisYear, ]


require(tidyr)

WDFWThisYearCounts <- WDFW_Crab_Logbooks_Data_ThisYearOnly
WDFWThisYearLbs <- WDFW_Crab_Logbooks_Data_ThisYearOnly
WDFWThisYearCounts<- WDFWThisYearCounts[ which(!is.na(WDFWThisYearCounts$Crab.Retained..count.)) ,]
WDFWThisYearLbs <- WDFWThisYearLbs[which(!is.na(WDFWThisYearLbs$Crab.Retained..lbs.)) ,]

WDFWThisYearCounts[,"CountsPerPot"] <- NA
WDFWThisYearLbs[,"CountsPerPot"] <- NA
WDFWThisYearCounts[,"LbsPerPot"] <- NA
WDFWThisYearLbs[,"LbsPerPot"] <- NA

#------------------

# CALCULATE COUNTS PER POT, AND LBS PER POT

WDFWThisYearCounts$CountsPerPot<-WDFWThisYearCounts$Crab.Retained..count./WDFWThisYearCounts$Pots.Fished
WDFWThisYearLbs$LbsPerPot<-WDFWThisYearLbs$Crab.Retained..lbs./WDFWThisYearLbs$Pots.Fished


#-------------
# cONSIDER THE COUNTS PER POT LOCATIONS
# But also convert to Pounds since that is what is used in Oregon and will be used for GAM.

LonsCounts<- -1*(WDFWThisYearCounts$Longitude.Begin.Degrees + WDFWThisYearCounts$Longitude.Begin.Minutes/60 )
LatsCounts<-WDFWThisYearCounts$Latitude.Begin.Degrees+ WDFWThisYearCounts$Latitude.Begin.Minutes/60

# CONVERT FROM COUNTS TO LBS FOR HAULS WHERE THIS IS NECESSARY. VERY ROUGH FOR NOW.
avgWtPerCrabLbs <- 1.8  # 1.8 from Carol Henry WDFW per comm May 4 2018

WDFWThisYearCounts$Crab.Retained..lbs.<-WDFWThisYearCounts$Crab.Retained..count.*avgWtPerCrabLbs

#--------------
# PLOT THE LBS PER POT LOCATIONS

LonsLbs<- -1*(WDFWThisYearLbs$Longitude.Begin.Degrees + WDFWThisYearLbs$Longitude.Begin.Minutes/60 )
LatsLbs<-WDFWThisYearLbs$Latitude.Begin.Degrees+ WDFWThisYearLbs$Latitude.Begin.Minutes/60


#--------------------
# NOW COMBINE THE LOGBOOK RECORDS that recorded things in lbs with the records that recorded things in counts.
# Note for now I assume no one recorded both lbs and counts in the same  logbook row
#-----------------

WDFWThisYearLbsAndCountsConvertedToLbs <-rbind(WDFWThisYearLbs, WDFWThisYearCounts)

WDFWThisYearLbsAndCountsConvertedToLbs$CountsPerPot<-WDFWThisYearLbsAndCountsConvertedToLbs$Crab.Retained..count./WDFWThisYearLbsAndCountsConvertedToLbs$Pots.Fished

#----------------------

Lons<- -1*(WDFWThisYearLbsAndCountsConvertedToLbs$Longitude.Begin.Degrees + WDFWThisYearLbsAndCountsConvertedToLbs$Longitude.Begin.Minutes/60 )
Lats<-WDFWThisYearLbsAndCountsConvertedToLbs$Latitude.Begin.Degrees+ WDFWThisYearLbsAndCountsConvertedToLbs$Latitude.Begin.Minutes/60

#--------------------
# PLOT THE POUNDS DATA THAT ULTIMATELY WE WILL USE FOR WASHINGTON
logbookPlotName<-paste(startYearForFishingSeason,"Washington Logbook Lbs")
PlotCoastForSardine(thisPlotName = logbookPlotName,myXlim=c(-126.0, -123.0), myYlim=c(45.0, 49.0))  # c(-126.0, -123.0),myYlim=c(42.5, 50.5)
title(main= logbookPlotName,cex.main=0.8)
points(Lons,Lats, pch=15,cex=0.2)
dev.copy(png,logbookPlotName)
dev.off()



#-----------------------

WDFWThisYearLbsAndCountsConvertedToLbs$Lats <- Lats
WDFWThisYearLbsAndCountsConvertedToLbs$Lons <- Lons

#----------------------
#----------------

require(geosphere)  # LOAD PACKAGE

#---------------

#----------------

for(i in 1:nrow(WDFWThisYearLbsAndCountsConvertedToLbs))
  {
 
  #-------------
   # For this logbook row's fishing month, subset oceanMatrix to only this month.
  calendarMonthNum <- WDFWThisYearLbsAndCountsConvertedToLbs$calendarMonthNum[i]
  rowsForThisMonth <- which(oceanMatrixThisYear[,"monthROMSOceanConditions"]==calendarMonthNum)

  # A quick check to make sure there's ROMS ocean data for this logbook point
   if(length(rowsForThisMonth)==0)
  {
    print(WDFWThisYearLbsAndCountsConvertedToLbs[i,])
    print("row")
    print(i)
    print("month")
     print(calendarMonthNum)
    stop("rows of ocean Jscope data for this fishing month are 0. Check if year (and months) match up in Jscope versus fishery data")
  }
  oceanMatrixThisYearThisMonth <- oceanMatrixThisYear[rowsForThisMonth, ]

  #----------
#  JAMEAL START HERE , AFTER LOADING THE GEOSPHERE PACKAGE. 

  ## reformat ocean lat,lon vectors into a matrix
oceanLonLatMatrix <- matrix(c(oceanMatrixThisYearThisMonth[,"nc_lonVector"], oceanMatrixThisYearThisMonth[ ,"nc_latVector"]), ncol = 2)

#  NOW FORMAT LON AND LAT
## subset lat,lon for single hake point
logbookLonLat <- c(WDFWThisYearLbsAndCountsConvertedToLbs$Lons[i],WDFWThisYearLbsAndCountsConvertedToLbs$Lats[i])



# THIS IS THE MAIN FUNCTION FROM MIKE MALICK TO GET DISTANCES BETWEEN LOGBOOK AND ROM POINTS
## calc dist between a single CRAB point and all ocean obs
dist.i <- geosphere::distHaversine(logbookLonLat, oceanLonLatMatrix)


## find index of ocean obs w/ min distance
ocean.ind <- which(dist.i == min(dist.i))[1]

## assign ocean obs w/ min distance to hake data

#-----------
# NOW ADD THE OCEAN VARIABLES AS NEW COLUMNS TO THE CRAB LOBGOOK DATA
#-----------
WDFWThisYearLbsAndCountsConvertedToLbs$distToNearestGridpointM[i] <- dist.i[ocean.ind]
WDFWThisYearLbsAndCountsConvertedToLbs$nc_lat[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "nc_latVector"]

#-----------
#   JUST PRINT SOME INFO TO SCREEN PERIODICALLY SO YOU KNOW IT IS RUNNING (SLOWLY)
if(i%%100==0)
{
print("row is")
print(i)
print("dist is")
print(WDFWThisYearLbsAndCountsConvertedToLbs$distToNearestGridpointM[i])
# just here for list of variables: oceanMatrix WAS cbind(nc_latVector,nc_lonVector,bathymetryVector,bottomTempVector,bottomOxyVector,bottomSaltVector,SSHVector,Chla2mVector,Chla10mVector,bottomAragVector,bottompHVector)
}
#---------------

}
#
#save("WDFWThisYearLbsAndCountsConvertedToLbs.Rdata")

return(WDFWThisYearLbsAndCountsConvertedToLbs)

}
