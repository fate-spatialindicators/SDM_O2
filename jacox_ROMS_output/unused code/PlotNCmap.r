PlotNCmap <- function (ncVariable, myTitle,nc_lon,nc_lat)
{
   # JUST PLOTS A GENERIC MAP OF US WEST COAST, USING OLD SARDINE MAPPING CODE 
  require(PBSmapping)
  source(here::here(
    "jacox_ROMS_output",
    "PlotCoastForSardine.r")
  )
  PlotCoastForSardine(myTitle)
  title(main= myTitle,cex.main=1.5)
  rbPal <-colorRampPalette(c('blue','red'))    # colorRampPalette(c('blue','red')
  
  # LOOK AT RANGE OF NCVARIABLE PASSED IN HERE
  thisMax<-max(ncVariable,na.rm = TRUE)
  thisMin<-min(ncVariable,na.rm = TRUE)
  thisMaxRounded <- round(thisMax,digits=1)
  thisMinRounded <- round(thisMin,digits=1)
  
   # COLOR IN THE MAP BASED ON VALUES FROM THE NCVARIABLE PASSED IN HERE. 
 thisRange<- (thisMax - thisMin)
  colorsForPredictions <- rbPal(10)[as.numeric(cut(ncVariable,breaks = seq(min(ncVariable,na.rm = TRUE),max(ncVariable,na.rm = TRUE),thisRange/10)))]
  points(nc_lon, nc_lat, pch=15,cex=0.2,col=colorsForPredictions)
  legend(x=-127.5,y=44.5,legend= c(thisMaxRounded,thisMinRounded),fill=c(rbPal(10)[10],rbPal(10)[1]))  # ,"Absent in survey","Predicted present","Predicted absent") ,fill=c(grayBlackPal(10)[10],grayBlackPal(10)[1],rbPal(12)[12],rbPal(12)[2]),bty="n")
  dev.off()
  
}