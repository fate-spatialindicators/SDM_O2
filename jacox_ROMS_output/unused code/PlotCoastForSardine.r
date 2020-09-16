PlotCoastForSardine <- function (thisPlotName,myXlim = c(-128.0, -123.0),myYlim=c(42.5, 50.5))

# Uses PBSmapping to plot polygons along the US West Coast, and color them in. 
# The idea is that below you read in EMOCCPolygonCoordsPBS, which defines (for example) 62 polygons along the US West Coast. 
# SEE THE PBSMAPPING DOCUMENTATION PDF FOR EXTENSIVE INFO ON THE PACKAGE AND AVAILABLE DATA !!!! 
#Isaac Kaplan 12/29/08

{
#  STOLEN FROM Code for Figure C1 of PBSmapping package documetation. 
#R: (Panel A)

#setwd('C:/Documents and Settings/kaplanis/My Documents/EconAtlantisProposal') # set working directory, place with  CSV file (see below)


require(PBSmapping)
palette("default")
data(nepacLL); # load the nepacLL data set
#data(worldLL)
png(filename=here::here(
  "jacox_ROMS_output",
  paste(thisPlotName,".png",sep="")
  )
)
plotMap(nepacLL, # plot the nepacLL data set
xlim=myXlim, # limit the region horizontally
ylim=myYlim, # limit the region vertically #plt=c(0.128, 0.776, 0.128, 0.776), # 0.16 0.97 OR 0.128 0.776 specify the plot region size
col=rgb(255, 255, 195, # set the foreground colour
maxColorValue=255),
bg=rgb(224, 253, 254, # set the background colour
maxColorValue=255),
tck=c(-0.03), # set the tick mark length
cex = 1.5, #was 1.8 adjust the font size
mgp=c(1.9, 0.4, 0),  # 1.9 0.7 0
border="black"); # adjust the axis label locations


# READ IN YOUR POLYGONS FROM CSV FILE
#EMOCCPolygonCoordsPBS<-read.table("EMOCCPolygonCoordsPBS.csv",as.is = TRUE,header=TRUE,sep=",")  # REad in data from Excel CSV file.  Headers are PID, POS, X, Y. Should be Polygon #, Point # within polygon, Long Dec Degress, Lat Dec degrees 
#EMOCCPolygonCoordsPBS<-as.PolySet(EMOCCPolygonCoordsPBS,projection="LL")  # Turn data into a "PolySet". The projection is "LL", i.e. lat-lon. 
#is.PolySet(EMOCCPolygonCoordsPBS) # just verify that the data is now a PolySet
}
