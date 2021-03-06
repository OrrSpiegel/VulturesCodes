#this code is genrated through  project from Rstudio on home PC newporject+GIT
#it will read an old mat file with hopefully gps data from the phd. 

# 22 Apr 2019
# Codes copied from LizardHRsizze2015.
#important params
HomePc=1; #for home PC 0 for lab
removeInterppoint=1;#to remove interpolated points?
SubsetData=0; #To set as 0 for not subsetting#the whole data  points is 985771 points for 2015. 

#D:\OrrS2\Box Sync\R codes\VulturesCodes\unitedGPSDataFromPhd.mat"

###### Helper functions ######
matlab2POS = function(x, timez = "UTC") {
  # Convert between MATLAB datenum values and R POSIXt time values.
  # 
  # Author: Luke Miller   Feb 20, 2011
  #Convert a numeric  MATLAB datenum (days since 0000-1-1 00:00) to seconds in 
  #the Unix epoch (seconds since 1970-1-1 00:00). Specify a time zone if the 
  #input datenum is anything other than the GMT/UTC time zone. 
  days = x - 719529   # 719529 = days from 1-1-0000 to 1-1-1970
  secs = days * 86400 # 86400 seconds in a day
  # This next string of functions is a complete disaster, but it works.
  # It tries to outsmart R by converting the secs value to a POSIXct value
  # in the UTC time zone, then converts that to a time/date string that 
  # should lose the time zone, and then it performs a second as.POSIXct()
  # conversion on the time/date string to get a POSIXct value in the user's 
  # specified timezone. Time zones are a goddamned nightmare.
  return(as.POSIXct(strftime(as.POSIXct(secs, origin = '1970-1-1', 
                                        tz = 'UTC'), format = '%Y-%m-%d %H:%M', 
                             tz = 'UTC', usetz = FALSE), tz = timez))
}

#plot.SpatialPolygons + downloaded helper functions- https://github.com/edzer/sp/blob/master/R/SpatialPolygons-displayMethods.R
{plot.SpatialPolygons <- function(x, col, border = par("fg"), add=FALSE, 
                                  xlim=NULL, ylim=NULL, xpd = NULL, density = NULL, angle = 45, 
                                  pbg=NULL, axes = FALSE, lty = par("lty"), ..., setParUsrBB=FALSE,
                                  usePolypath=NULL, rule=NULL, bgMap = NULL) {
  
  if (is.null(pbg))
    pbg = par("bg") # transparent!
  if (!is(x, "SpatialPolygons")) 
    stop("Not a SpatialPolygons object")
  if (is.null(usePolypath)) usePolypath <- get_Polypath()
  if (is.null(rule)) rule <- get_PolypathRule()
  
  if (! add) 
    plot(as(x, "Spatial"), xlim=xlim, ylim=ylim, axes = axes, 
         ..., setParUsrBB=setParUsrBB, bgMap = bgMap)
  
  n <- length(slot(x, "polygons"))
  if (length(border) != n)
    border <- rep(border, n, n)
  polys <- slot(x, "polygons")
  pO <- slot(x, "plotOrder")
  if (!is.null(density)) {
    if (missing(col)) col <- par("fg")
    if (length(col) != n) col <- rep(col, n, n)
    if (length(density) != n)
      density <- rep(density, n, n)
    if (length(angle) != n)
      angle <- rep(angle, n, n)
    for (j in pO) 
      .polygonRingHoles(polys[[j]], border = border[j], 
                        xpd = xpd, density = density[j], angle = angle[j], 
                        col = col[j], pbg = pbg, lty=lty, ...) 
  } else {
    if (missing(col)) col <- NA
    if (length(col) != n) col <- rep(col, n, n)
    for (j in pO) 
      .polygonRingHoles(polys[[j]], col=col[j], 
                        border=border[j], xpd = xpd, pbg = pbg, lty=lty, ...,
                        usePolypath=usePolypath, rule=rule)
  }
}
  
  setMethod("plot", signature(x = "SpatialPolygons", y = "missing"),
            function(x, y, ...) plot.SpatialPolygons(x, ...))
  
  .polygonRingHoles <- function(Sr, col=NA, border=NULL, xpd=NULL, density=NULL,
                                angle=45, pbg, lty = par("lty"), ..., usePolypath=NULL,
                                rule=NULL) {
    if (!is(Sr, "Polygons")) 
      stop("Not an Polygons object")
    if (is.null(usePolypath)) usePolypath <- get_Polypath()
    if (is.null(rule)) rule <- get_PolypathRule()
    if (!is.null(density)) hatch <- TRUE
    else hatch <- FALSE
    pO <- slot(Sr, "plotOrder")
    polys <- slot(Sr, "Polygons")
    
    if (hatch) {
      for (i in pO) {
        if (!slot(polys[[i]], "hole"))
          .polygon(slot(polys[[i]], "coords"), 
                   border = border, xpd = xpd, 
                   density = density, angle = angle,
                   col=col, hatch=TRUE, lty=lty, ...)
        else .polygon(slot(polys[[i]], "coords"), 
                      border = border, xpd = xpd, col=pbg, 
                      density = NULL, lty=lty, ...)
      } 
    } else if (exists("polypath") && usePolypath) {
      Srl <- as(Sr, "Lines")
      crds <- coordinates(Srl)
      if (length(crds) == 1) mcrds <- crds[[1]]
      else {
        NAr <- as.double(c(NA, NA))
        crds1 <- lapply(crds, function(x) rbind(x, NAr))
        mcrds <- do.call(rbind, crds1)
        mcrds <- mcrds[-nrow(mcrds),]
        rownames(mcrds) <- NULL
      }
      polypath(x=mcrds[,1], y=mcrds[,2], border=border, col=col,
               lty=lty, rule=rule, xpd=xpd, ...)
    } else {
      for (i in pO) {
        if (!slot(polys[[i]], "hole"))
          .polygon(slot(polys[[i]], "coords"), 
                   border = border, xpd = xpd, 
                   col=col, lty=lty, ...)
        else .polygon(slot(polys[[i]], "coords"), 
                      border = border, xpd = xpd, col=pbg, lty=lty,
                      ...)
      }
    }
  }
  
  
  .polygon = function(x, y = NULL, density = NULL, angle = 45,
                      border = NULL, col = NA, lty = NULL, xpd = NULL, hatch=NA, ...) {
    if (is.na(hatch)) polygon(x = x, y = y, border = border, 
                              col = col, lty = lty, xpd = xpd, ...)
    else polygon(x = x, y = y, density = density, angle = angle, 
                 border = border, lty = lty, xpd = xpd, col=col, ...)
  }
}



#Required packages #######
require(R.matlab);require(CircStats); #require(boot) ; require(MASS) ;   require(fields); 
require(adehabitatHR);require(plyr);#require(lme4);require(AICcmodavg);require(ggplot2);
require(rgeos);require(sp);#require(spatstat);   require(scales);



##### Loading Vultures Data from an old mat file and some filtering ####
if (HomePc==1) {
  m.input.path <- "D:\\OrrS2\\Box Sync\\R codes\\VulturesCodes\\" 
  # now as a part of a project setwd("D:\\Dropbox\\R codes\\HRandMovement codes")
  #r.out.path <- "D:\\Dropbox\\Matlab\\Lizards\\matFiles\\fromR\\" 
}else if (HomePc==0) {
  m.input.path <- "C:\\Users\\Ors\\OneDrive\\MATLAB\\Lizards\\MatFiles\\" 
  #setwd("C:\\Users\\ors\\Dropbox\\R codes\\HRandMovement codes")
  #r.out.path <- "D:\\Dropbox\\Matlab\\Lizards\\matFiles\\fromR\\" 
}


#files = list.files(m.input.path, full.names=TRUE)
#if(! (length(files)>0)){print("no files in folder or wrong path")}  
File.name="UnitedGPSDataFromPhD"; #
Data.FromMat <-readMat(paste(m.input.path,File.name,".mat", sep = ""))
Data.FromMat[['GeoStruct']]=NULL
View(Data.FromMat[['GeoStruct']])
dataStrct=Data.FromMat[["TagsDataStrctrIntrpl"]]
save(dataStrct,file='temp.rdata')

str(Data.FromMat)
#columns in 1=tag 2=date&time 3=date only 4=time only 5=Lat 6=Long 7=elevation 8=speed 
#        9=Numof Sattlt 10=DOP param 11=TGSV_val 12=UTM_Easting 13=UTM_Northing 
#        14 distance from center 15 burst number, 16 PointTimeDiff 17  PointDistDiff
#        18 CalcSpeed %19 interpolated


##### processing mat file and saving it ######
#converting into a data.frame and giving names, changing date format
MetaDataOnLiz=data.frame(Data.FromMat$MetaDataOnLiz)
names(MetaDataOnLiz)=c('PropOfFilteredFix','DaysWithData','TotalBursts','AvgBurstsperday','propfilteredpoints1','filteredpoints2','firstDay','LastDay')
MetaDataOnLiz$firstDay=matlab2POS(MetaDataOnLiz$firstDay);#converting dates to r format
MetaDataOnLiz$LastDay=matlab2POS(MetaDataOnLiz$LastDay);#converting dates to r format

LizMatData=data.frame(Data.FromMat$MergedLizMatrixByBurst)
LizMatData[LizMatData=="NaN"]<-NA#replacing NaN with NA, takes time.
names(LizMatData)=c('names','dateWtime','date_only','time_only','Lat','Long','elevation','speed','NumofSattlt',
                    'DOP_param','TGSV_val','UTM_Easting','UTM_Northing','DistFromCenter','burst_number','PointTimeDiff',
                    'PointDistDiff','CalcSpeed','interpolated')

#saving days as burst with a different name for each indiv
LizMatData$DaysAsBurst=as.factor(paste(1+LizMatData$date_only-min(LizMatData$date_only),LizMatData$names,2015,sep="_"));

LizMatData$POSIXctTime=matlab2POS(LizMatData$dateWtime);#converting dates to r format
LizMatData$date_only=matlab2POS(LizMatData$date_only);#converting dates to r format
LizMatData$time_only=matlab2POS(1+LizMatData$time_only);#converting dates to r format must be more than 1 for min date

LizMatData$NameFactor= as.factor(paste(LizMatData$names,2015,sep="_"))#name_year as factor
LizMatData$names= as.factor(LizMatData$names)#name as factor

#throwing all interpolated points  if removeInterppoint= 1 deleting interpolated points #
if (removeInterppoint==1){LizMatData=LizMatData[ which(LizMatData$interpolated==0), ]}#for 2015 no interpolation was done anyway

head(LizMatData)
rm(Data.FromMat, m.input.path)

## saving the whole dataset
Name='Liz_2015GPS_Data_R_FormatAll.rdata'
save(list=c('LizMatData','MetaDataOnLiz'),file=Name)


#with(Data.FromMat,data.frame(names,years,Lines,GPSlinesAfterFilter,MissingDaysPerLizard, DaysRange));
#LizMatData=with(Data.FromMat,data.frame(UUTM.Easting ,UUTM.Northing,StepUnHorizAcc,StepUHorizDilPrec,UCalcSpeed,UGPSDateAsNum,UGPSTimeAsNum,UGPSyear,UGPSname,UPointDistDiff,UPointTimeDiff,UisInterpolated,Lat,Lon,UDayFirstPoint,UDayLastPoint,DistFromRoad,DistFromCenter,MindistFromDams,UObsthreatCont,UCurisotyCont,UConspesAgress,RecentSteps));


#which(tt$StepUnHorizAcc==NaN)
#str(LizMatData); View(LizMatData)


# ##### working on a substet for the HR analsis set  SubsetData to 0 to skip #####  
#LizMatData2=LizMatData[seq(from=1, to=100000,by=100),]
if (SubsetData>0){#working on data subset for the HR analysis
  LizMatData2=LizMatData[seq(from=1,to=dim(LizMatData)[1],by=SubsetData),c(1,12,13,20,21,22)]#the subset has also only some of the columns!!!
  head(LizMatData2)
  LizMatData2 <- droplevels(LizMatData2); 
  #   #name_year as factor
  #   LizMatData2$NameFactor= as.factor(paste(LizMatData2$UGPSname,LizMatData2$UGPSyear,sep="_"))
  #   #saving days as burst with a different name for each indiv
  #   LizMatData2$DaysAsBurst=as.factor(paste(1+LizMatData2$UGPSDateAsNum-min(LizMatData2$UGPSDateAsNum),LizMatData2$UGPSname,LizMatData2$UGPSyear,sep="_"));
  #   
}else{LizMatData2=LizMatData}#working on a substet


##### as.ltraj class #####
print('now converting to track ltraj ')
#conversion to class of adehabitatLT
LizXYdataSubset=as.ltraj(xy = LizMatData2[,c("UTM_Easting","UTM_Northing")], 
                         date = LizMatData2$POSIXctTime, 
                         id = LizMatData2$NameFactor,
                         burst = LizMatData2$DaysAsBurst,
                         typeII=TRUE,
                         infolocs =LizMatData2);

LizXYdataFull=as.ltraj(xy = LizMatData[,c("UTM_Easting","UTM_Northing")], 
                       date = LizMatData$POSIXctTime, 
                       id = LizMatData$NameFactor,
                       burst = LizMatData$DaysAsBurst,
                       typeII=TRUE,
                       infolocs =LizMatData);
#in addiotn the infolocs include the distance and R2N and dx dy and rel.angle turnign angle
head(LizXYdataSubset[[2]])
#just to work on a small file with bursts (days) from 10 lizards
LizXYdataSample=LizXYdataSubset[id=id(LizXYdataSubset[10])]
#LizXYdataSample=LizXYdataFull  [id=id(LizXYdataFull)[1:10]]

#this will plot the tracks of all lizards. heavy
#plot(LizXYdata)
#plot(LizXYdataSample)

## Saving processed file first time, before HR analysis ####
Name='Liz_2015GPS_Data_R_FormatW_xy.rdata'
save(list=ls(),file=Name)
print('saved  Lizard GPS data  ')





