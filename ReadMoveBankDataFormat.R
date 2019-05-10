# Apr 23rd :-) 2019 trying to read movebank data using move


#### getting packages #######
require(move)#for downloading data
require(mapproj);require(ggmap) #these packages are necessary to work with google maps
require(spatsoc);require("asnipe");require("igraph"); # for working with the SN parts

#require(dplyr) 
#require(tidyverse) 

#### key paramterer values ######
MaxSpeedPermited=120 #in movebank units (m/s??) anything above this will be filtered
MB.LoginObject=movebankLogin(username='ors',password='eesH1eeF')
VulturesToPlotMap=10:15 #1:length(unstackedOhad) #out of the vultures in the DB which one to plot? choose a few out of the 83
DistThresholM=2000 #in meters
TimeThreshold='10 minutes' #in this format 'XX units'

################### reading data from movebank ###################
#movebank user ors PW: eesH1eeF
#"Gyps fulvus INPA Hatzofe" Movebank ID	6071688
#"HUJ MoveEcol Lab Israel: Griffon vulture Gyps fulvus" Movebank ID	6638215
#"E-obs Gsm Vultures Israel" Movebank ID	7359070
#timestamp_start, timestamp_end character or POSIXct=?yyyyMMddHHmmssSSS?
MoveStackDatasetOhad=getMovebankData(study=6071688, login=MB.LoginObject,
        includeExtraSensors=FALSE, deploymentAsIndividuals=FALSE,removeDuplicatedTimestamps=TRUE)#animalName
MoveStackDatasetOrr=getMovebankData(study=6638215, login=MB.LoginObject,
                                 includeExtraSensors=FALSE, deploymentAsIndividuals=FALSE,removeDuplicatedTimestamps=TRUE)#animalName

#returns or a MoveStack object can use as.data.frame
#corridor
unstackedOhad=split(MoveStackDatasetOhad) #splitting Ohads data into individuals 
#summary(unstackedOhad)  


#### just loooking on raw downloaded data ####
show(MoveStackDatasetOhad)
summary(MoveStackDatasetOrr)
citations(MoveStackDatasetOrr)
equalProj(MoveStackDatasetOrr)
n.locs(MoveStackDatasetOhad)
timeLag(MoveStackDatasetOhad)
#too heavy but works: plot(MoveStackDatasetOhad,  lwd=2, xlab="location_long", ylab="location_lat")


############# very basic filtering. stage 1  ##############
#without going to a data frame since i cannot go back and doing it one by one (another try)###
## creating an empty metadata storage df
TagsMetaData=setNames(data.frame(matrix(ncol =6, nrow = length(unstackedOhad))), c('name',"InitialPoints", "PercentThrown", "N_locs","TrackDurationDays","TrackDurationStartEnd"))

## a loop on tags. Need to check tags: 75,76,77,78, and maybe also 3 and 10. ###### #
for (indv in 1:length(unstackedOhad) ){## loop on individuals, now in separate Move objects
  TagsMetaData$InitialPoints[indv]=  dim(unstackedOhad[[indv]]@data)[1]
  TagsMetaData$name[indv]=          names(unstackedOhad)[indv]
  plot(unstackedOhad[[indv]],col='blue', type='b',main=paste("Indiv=",indv, ', name=',TagsMetaData$name[indv],sep=' '))#just simple plot of this individual's points
  #dim(unstackedOhad[[indv]]@coords)
  
  ## removing unneeded columns
  VarsToRemove <- names(unstackedOhad[[indv]]@data) %in% c("sensor_type_id","taxon_canonical_name","nick_name","earliest_date_born","sensor","optional",
                                            "sensor_type","mw_activity_count","eobs_accelerations_raw","eobs_acceleration_sampling_frequency_per_axis",
                                            "eobs_acceleration_axes","argos_valid_location_algorithm","argos_sensor_4","argos_sensor_3","argos_sensor_2",
                                            "argos_sensor_1","argos_semi_minor","argos_semi_major","argos_pass_duration","argos_orientation","argos_nopc",
                                            "argos_lat1","argos_lat2","1084088","argos_lon1","argos_lon2","argos_nb_mes","argos_nb_mes_120",
                                            "eobs_key_bin_checksum","eobs_fix_battery_voltage","eobs_battery_voltage","eobs_status",
                                            "eobs_start_timestamp","eobs_type_of_fix","eobs_used_time_to_get_fix","eobs_temperature",
                                            "gps_dop","magnetic_field_raw_x","magnetic_field_raw_y","magnetic_field_raw_z","ornitela_transmission_protocol",
                                            "tag_voltage","algorithm_marked_outlier","argos_altitude","argos_best_level","argos_lc","argos_iq",
                                            "argos_gdop","argos_error_radius","argos_calcul_freq","location_lat.1","location_long.1","timestamps","height_raw",
                                            "barometric_pressure","barometric_height","battery_charging_current","eobs_activity","manually_marked_outlier",
                                            "eobs_activity_samples", "acceleration_raw_y", "battery_charge_percent", "data_decoding_software","gps_vdop","height_above_ellipsoid",
                                            'acceleration_raw_x','acceleration_raw_z',"acceleration_raw_z","eobs_horizontal_accuracy_estimate","eobs_speed_accuracy_estimate");  
  unstackedOhad[[indv]]@data=unstackedOhad[[indv]]@data[!VarsToRemove]
  dim(unstackedOhad[[indv]]@data)#colunms removed?
  
  ## filtering: choosing indices to keep for this individual
  
  indx=1:dim(unstackedOhad[[indv]]@data)[1] #starting with all points a
  if(sum(unstackedOhad[[indv]]@data$heading <= 360,na.rm=T)){#do i have heading data or this one?
    indx=intersect(indx,which(unstackedOhad[[indv]]@data$heading <= 360))} #if yes, now index include only points with realistic heading 
  if(sum(unstackedOhad[[indv]]@data$ground_speed<=MaxSpeedPermited,na.rm=T)){#below threshhold speed?
    indx=intersect(indx,which(unstackedOhad[[indv]]@data$ground_speed<=120))}
  #if(sum(unstackedOhad[[indv]]@data$gps_time_to_fix<=89,na.rm=T)){
  #   indx=intersect(indx,which(unstackedOhad[[indv]]@data$gps_time_to_fix<=89))}
  if(sum(unstackedOhad[[indv]]@data$gps_satellite_count>=3,na.rm=T)){#enough satellite numbers?
    indx=intersect(indx,which(unstackedOhad[[indv]]@data$gps_satellite_count>=3))}
  
  
  ## subsetting the different slots of this move object
  print(paste("indiv",indv,"name",TagsMetaData$name[indv],'. I throw out', TagsMetaData$InitialPoints[indv]-length(indx), 'points, out of',TagsMetaData$InitialPoints[indv]))
  TagsMetaData$PercentThrown[indv]=(TagsMetaData$InitialPoints[indv]-length(indx))/TagsMetaData$InitialPoints[indv]
  
  unstackedOhad[[indv]]@timestamps=unstackedOhad[[indv]]@timestamps[indx]
  #unstackedOhad[[indv]]@idData
  unstackedOhad[[indv]]@sensor=unstackedOhad[[indv]]@sensor[indx]
  unstackedOhad[[indv]]@data=unstackedOhad[[indv]]@data[indx,]
  #unstackedOhad[[indv]]@coords.nrs
  unstackedOhad[[indv]]@coords=unstackedOhad[[indv]]@coords[indx,]
  unstackedOhad[[indv]]@bbox[1,]=range(unstackedOhad[[indv]]@coords[,1]);unstackedOhad[[indv]]@bbox[2,]=range(unstackedOhad[[indv]]@coords[,2])
  #unstackedOhad[[indv]]@proj4string
  
  
  ## collecting metadata and plotting fitered track: 
  TagsMetaData$N_locs[indv]=  dim(unstackedOhad[[indv]]@data)[1]
  TagsMetaData$TrackDurationDays[indv]=  length(unique(as.Date(as.character(unstackedOhad[[indv]]@data$timestamp))))
  TagsMetaData$TrackDurationStartEnd[indv]=  (max(as.Date(as.character(unstackedOhad[[indv]]@data$timestamp)))-min(as.Date(as.character(unstackedOhad[[indv]]@data$timestamp))))
  lines(unstackedOhad[[indv]],col='red')
  #plot(unstackedOhad[[indv]], type="o", col=3, lwd=2, pch=20, xlab="location_long", ylab="location_lat")
  
  ##logging metadata
  
  head(timeLag(unstackedOhad[[indv]], units="mins"))
  head(timestamps(unstackedOhad[[indv]]))
  
}#loop on individuals


## merging back to a unifide movestack object
BackStackedOhad=moveStack(unstackedOhad) #merging for next section Ohads data into individuals 


#### mapping all the filtered dataset together: #####
#preparing a dataframe for mapping:
DatasetOhadF=as.data.frame(BackStackedOhad) #converting to a dataset 
DatasetOhadF$DateOnly=as.Date(as.character(DatasetOhadF$timestamp));

mapObj <- get_map(bbox(extent(min(DatasetOhadF$location_long), 
                              max(DatasetOhadF$location_long), 
                              min(DatasetOhadF$location_lat), 
                              max(DatasetOhadF$location_lat) )*1.1), source="google", zoom=6)
## plot 1 colors by time
Map=ggmap(mapObj)+
  geom_path(data=DatasetOhadF, aes(x=location_long, y=location_lat,colour = as.numeric(DateOnly)), size=1) +
  scale_colour_gradient2(low = "blue",mid='purple',high = "red",midpoint = mean(as.numeric(DatasetOhadF$DateOnly)))+
  #ggtitle(paste("Indiv=",indv, ', name=',TagsMetaData$name[indv],"duration",TagsMetaData$TrackDurationDays[indv],'day',sep=' '))
  ggtitle('all tags together, color as time gradient');
print(Map);rm(Map)

## plot 2 colors by tag
Map=ggmap(mapObj)+
  geom_path(data=DatasetOhadF, aes(x=location_long, y=location_lat,colour = ID), size=1,show.legend = FALSE) +
  #scale_colour_gradient2(low = "blue",mid='purple',high = "red",midpoint = mean(as.numeric(DatasetOhadF$DateOnly)))+
  #ggtitle(paste("Indiv=",indv, ', name=',TagsMetaData$name[indv],"duration",TagsMetaData$TrackDurationDays[indv],'day',sep=' '))
  ggtitle('all tags together, colors but ID');print(Map)
print(Map);rm(Map)


### mapping by tags  #####
## another loop on indiv for plotting a background map
for (indv in VulturesToPlotMap ){## loop on individuals, now for plotting
  ## converting to a df and adding a date only
  DFAsingleIndiv=as.data.frame(unstackedOhad[[indv]])
  DFAsingleIndiv$DateOnly=(as.Date(as.character(DFAsingleIndiv$timestamp)));
   
  #mapObj <- get_map(bbox(extent(DFAsingleIndiv)*1.1), source="stamen", zoom=12) #something is still fucked with the extent 
  # mapObj <- get_map(bbox(extent(min(DFAsingleIndiv$location_long), 
  #                          max(DFAsingleIndiv$location_long), 
  #                          min(DFAsingleIndiv$location_lat), 
  #                          max(DFAsingleIndiv$location_lat) )*1.1), source="stamen", zoom=10)
  
  mapObj <- get_map(bbox(extent(min(DFAsingleIndiv$location_long), 
                                max(DFAsingleIndiv$location_long), 
                                min(DFAsingleIndiv$location_lat), 
                                max(DFAsingleIndiv$location_lat) )*1.1), source="google", zoom=8)
  #plots with a fixed color
  # Map=ggmap(mapObj)+
  #   geom_path(data=DFAsingleIndiv, aes(x=location_long, y=location_lat),colour="red", size=1) +
  #   ggtitle(paste("Indiv=",indv, ', name=',TagsMetaData$name[indv],sep=' '))
  # print(Map)
  # rm(Map)
  # 
  Map=ggmap(mapObj)+
    geom_path(data=DFAsingleIndiv, aes(x=location_long, y=location_lat,colour = as.numeric(DateOnly)), size=1) +
    scale_colour_gradient2(low = "blue",mid='purple',high = "red",midpoint = mean(as.numeric(DFAsingleIndiv$DateOnly)))+
    ggtitle(paste("Indiv=",indv, ', name=',TagsMetaData$name[indv],"duration",TagsMetaData$TrackDurationDays[indv],'day',sep=' '))
  
  print(Map)
  rm(Map)
  
  }#loop

####################### towards an adjacency matrix #############
## converting to a df and than data.table after fitlering
#now done above: BackStackedOhad=moveStack(unstackedOhad) #splitting Ohads data into individuals 
#too heavy but works: plot(BackStackedOhad,  lwd=2, xlab="location_long", ylab="location_lat")
DatasetOhadF=as.data.frame(BackStackedOhad) #converting to a dataset 
DatasetOhadF=data.table::setDT(DatasetOhadF,keep.rownames=T)#and now to a data.table

## renaming\adding columns
DatasetOhadF$ID=DatasetOhadF$trackId
DatasetOhadF$location_long2=DatasetOhadF$location_long
DatasetOhadF$location_lat1=DatasetOhadF$location_lat

## adding coordinates with a metric value ##### #
## first setting the same coordinate system as in movebank
projection(MoveStackDatasetOhad)#this was the original data projection from movebank
DatasetOhadF_wgs=DatasetOhadF;
coordinates(DatasetOhadF_wgs)<-~location_long+location_lat
proj4string(DatasetOhadF_wgs)<-CRS(projection(MoveStackDatasetOhad))

#proj4string(x) <-CRS("+proj=utm +zone=10+datum=WGS84")
#newData <- spTransform(x, CRS("+init=epsg:4238"))
utmS <- '+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs' #south, but most points are in the north
utmN <- '+proj=utm +zone=36        +ellps=WGS84 +datum=WGS84 +units=m +no_defs'  #north

## converting to UTM36North, but note that not all points are within, just the majority
DatasetOhadF_utmN <- spTransform(DatasetOhadF_wgs, CRS(utmN))
head(coordinates(DatasetOhadF_utmN))#now the lat long are in metric
# just plotting the new dataset to see the locations look fine: plot(DatasetOhadF_utmN,col='blue')

## appending the new coordinates in easting northing
DatasetOhadF$Easting=coordinates(DatasetOhadF_utmN)[,1]
DatasetOhadF$Northing=coordinates(DatasetOhadF_utmN)[,2]

## groupping into time groups
group_times(DatasetOhadF, datetime = 'timestamp', threshold = TimeThreshold)

## grouping into spatial groups with the metric values
group_pts(DatasetOhadF, threshold = DistThresholM, id = 'ID', coords = c('Northing', 'Easting'), timegroup = 'timegroup')
#use group_lines instead??

#group_times(DT = DT, datetime = 'datetime', threshold = '1 day')
#group_lines(DT, threshold = 50, projection = utm,   id = 'ID', coords = c('X', 'Y'),timegroup = 'timegroup', sortBy = 'datetime')
edge_nn(DatasetOhadF, id = 'ID', coords = c('Northing', 'Easting'), timegroup = 'timegroup')
View(DatasetOhadF)

#### 
gbiOhad <- get_gbi(DatasetOhadF, group = 'group', id = 'ID')
netOhad <- get_network(gbiOhad, data_format = "GBI", association_index = "SRI")

g <- graph.adjacency(netOhad, 'undirected', diag = FALSE, weighted = TRUE)

##### manual distanc calculation and SN construction#####
CoocuranceCounter   <- expand.grid(unique(as.character(DatasetOhadF$ID)),unique(as.character(DatasetOhadF$ID)))#a long form of all possible dyads to count interaction
SimDatapointCounter <- expand.grid(unique(as.character(DatasetOhadF$ID)),unique(as.character(DatasetOhadF$ID)))#a long form of all possible dyads to count intervals both were atthe same timegroup


ColumToSelect=c("ID","Northing","Easting","location_lat","location_long",'timegroup',"group")
for (timgrpind in 1: max(DatasetOhadF$timegroup)){
  ## extract current time group #18458 has a good example
  #subset(DatasetOhadF, timegroup==timgrpind,select=c("ID","location_lat","location_long"))
  timegroupDF=subset(DatasetOhadF, timegroup==timgrpind,select=ColumToSelect)
  
  if (dim(timegroupDF)[1]>1){#more than one vulture in this time group?
    print(timegroupDF[order(timegroupDF$group),])#,'timegroupDF$ID'#plot it
    
    # create all possible dyads of vultures in a long format
    timegroupDF=subset(timegroupDF, timegroup==timgrpind,select=c("ID","location_lat","location_long","group"))
    dt <- expand.grid.df(timegroupDF,timegroupDF)
    names(dt)[6:7] <- c("lat_secondInd","long_secondInd")
    
    setDT(dt)[ , dist_km := distGeo(matrix(c(location_long, location_lat), ncol = 2), 
                                    matrix(c(long_secondInd, lat_secondInd), ncol = 2))/1000];
    
    ##subsetting to the witihn distancethreshold lines
    PresentVultures=subset(dt, dist_km==0,select=c(1,5))#these are the vultures present in this time group
    as.character(PresentVultures[,1])==as.character(PresentVultures[,2])
    #PresentVultures=unique(PresentVultures$ID);
    
    InteractingDyads=subset(dt, dist_km<=(DistThresholM/1000)&dist_km!=0)
    print(InteractingDyads[order(InteractingDyads$group,InteractingDyads$ID),])#,'timegroupDF$ID'#plot it
    
    rm(dt)
  }#more than one vulture in this time group?
  
}#loop on time groups



#### Calculation of the utilization distribution ####
brownian.bridge.dyn




########## DRAFTS #############################################################################################################

#### conversion from Move Class to DF and removal of non needed columns ####
DatasetOrr=as.data.frame(MoveStackDatasetOrr) #converting to a dataset 
DatasetOhad=as.data.frame(MoveStackDatasetOhad) #converting to a dataset 
#str(DatasetOrr) ;str(DatasetOhad)
## removing un-necessary columns
unique(DatasetOrr$deployment_id)
VarsToRemove <- names(DatasetOrr) %in% c("sensor_type_id","update_ts","sensor_type","sensor","earliest_date_born",
                                         "exact_date_of_birth", "latest_date_born","nick_name","ring_id","taxon_canonical_name",
                                         "location_lat.1","location_long.1","timestamps","comments.1","death_comments");DatasetOrr=DatasetOrr[!VarsToRemove]
View(DatasetOrr)
#"argos_best_level","argos_error_radius","argos_iq",""
DatasetOhad=subset(DatasetOhad, sensor_type %in% c("GPS") )#geeting rid of argos lines
unique(DatasetOhad$sensor_type_id)
VarsToRemove <- names(DatasetOhad) %in% c("sensor_type_id","taxon_canonical_name","nick_name","earliest_date_born","sensor","optional",
                                          "sensor_type","mw_activity_count","eobs_accelerations_raw","eobs_acceleration_sampling_frequency_per_axis",
                                          "eobs_acceleration_axes","argos_valid_location_algorithm","argos_sensor_4","argos_sensor_3","argos_sensor_2",
                                          "argos_sensor_1","argos_semi_minor","argos_semi_major","argos_pass_duration","argos_orientation","argos_nopc",
                                          "argos_lat1","argos_lat2","1084088","argos_lon1","argos_lon2","argos_nb_mes","argos_nb_mes_120",
                                          "eobs_key_bin_checksum","eobs_fix_battery_voltage","eobs_battery_voltage","eobs_status",
                                          "eobs_start_timestamp","eobs_type_of_fix","eobs_used_time_to_get_fix","eobs_temperature",
                                          "gps_dop","magnetic_field_raw_x","magnetic_field_raw_y","magnetic_field_raw_z","ornitela_transmission_protocol",
                                          "tag_voltage","algorithm_marked_outlier","argos_altitude","argos_best_level","argos_lc","argos_iq",
                                          "argos_gdop","argos_error_radius","argos_calcul_freq","location_lat.1","location_long.1","timestamps","height_raw",
                                          "barometric_pressure","barometric_height","battery_charging_current","eobs_activity","manually_marked_outlier",
                                          "eobs_activity_samples", "acceleration_raw_y", "battery_charge_percent", "data_decoding_software","gps_vdop","height_above_ellipsoid",
                                          'acceleration_raw_x','acceleration_raw_z',"acceleration_raw_z","eobs_horizontal_accuracy_estimate","eobs_speed_accuracy_estimate");  
DatasetOhad=DatasetOhad[!VarsToRemove]


#DatasetOrr$update_ts=NULL;DatasetOrr$taxon_canonical_name=DatasetOrr$TBA=NULL #Orr

#myvars <- names(mydata) %in% c("v1", "v2", "v3")
#unique(DatasetOrr$sensor_type_id)
unique(DatasetOhad$gps_satellite_count)
View(DatasetOhad)
#DatasetOhadF$height_above_ellipsoid

### initial dataset exploration & simple filteres, but this has to be modified to do it one by one since not all have all data colums ######
range(DatasetOhad$ground_speed,na.rm=T)#needs filtering

DatasetOhadF <- subset(DatasetOhad, heading <= 360 & gps_time_to_fix<=89 & ground_speed<=120 & gps_satellite_count>=3  ) #,select=c(ID, Weight))
#hist(DatasetOhadF$gps_time_to_fix,1000);
#hist(DatasetOhadF$gps_time_to_fix[DatasetOhadF$gps_time_to_fix>70&DatasetOhadF$gps_time_to_fix<100],100)
#hist(DatasetOhadF$ground_speed[DatasetOhadF$ground_speed>40],100)
range(DatasetOhadF$gps_time_to_fix,na.rm=T) #now ok, threw out everything above 89 since there is a peack in 90
unique(DatasetOhadF$gps_satellite_count)
range(DatasetOhadF$height_above_ellipsoid,na.rm=T) #now ok, threw out everything above 89 since there is a peack in 90
hist(DatasetOhadF$location_lat,100)
n.locs(DatasetOhadF)





#a code from 
library(sf)
library(tidyverse)

# Coordinate examples with expected UTM values
coord_sample <- data.frame(
  "Northing" = c(1105578.589, 5540547.370),
  "Easting" = c(609600.773, 643329.124),
  "Latitude" = c(10, 50),
  "Longitude" = c(118, 119))

#' Find UTM EPSG code from Latitude and Longitude coordinates (EPSG 4326 WGS84)
#from trunning the code i find that EPSG
#code:   32633   32634   32635   32636   32637 
#n_locs: 3769     847     750   1071581   15594 
LatLonToUTMEPSGCode <- function(lat, lon) {
#' (vectorised)
#' Source: https://geocompr.robinlovelace.net/reproj-geo-data.html
#' Source: https://gis.stackexchange.com/questions/13291/computing-utm-zone-from-lat-long-point
  
  zone_number <- (floor((lon + 180) / 6) %% 60) + 1
  
  # Special zones for Norway
  cond_32 <- lat >= 56.0 & lat < 64.0 & lon >= 3.0 & lon < 12.0
  zone_number[cond_32] <- 32
  
  # Special zones for Svalbard
  cond_lat <- lat >= 72.0 & lat < 84.0
  
  cond_31 <- cond_lat & lon >= 0.0 & lon <  9.0
  zone_number[cond_31] <- 31
  
  cond_33 <- cond_lat & lon >= 9.0 & lon < 21.0
  zone_number[cond_33] <- 33
  
  cond_35 <- cond_lat & lon >= 21.0 & lon < 33.0
  zone_number[cond_35] <- 35
  
  cond_37 <- cond_lat & lon >= 33.0 & lon < 42.0
  zone_number[cond_37] <- 37
  
  # EPSG code
  utm <- zone_number
  utm[lat > 0] <- utm[lat > 0] + 32600
  utm[lat <= 0] <- utm[lat <= 0] + 32700
  
  return(utm)
}

sf_sample <- sf::st_as_sf(coord_sample , coords = c("Longitude", "Latitude"), crs = 4326)
sf_sample2 <- sf::st_as_sf(DatasetOhadF, coords = c("location_long", "location_lat"),crs = 4326)
#EPSG=32636  
#proj4 +proj=utm +zone=36 +ellps=WGS84 +datum=WGS84 +units=m +no_defs 

sf_sample %>%
  do(cbind(., st_coordinates(.))) %>%
  mutate(EPSG = LatLonToUTMEPSGCode(lat = Y, lon = X)) %>%
  group_by(EPSG) %>%
  do(cbind(as.data.frame(.) %>% select(Northing, Easting),
           st_coordinates(st_transform(., crs = head(.$EPSG, 1))))) %>%
  ungroup()

EPSG = LatLonToUTMEPSGCode(lat =DatasetOhadF$location_lat, lon = DatasetOhadF$location_lat)





library(maps)
library(reshape)  
library(dplyr)
library(data.table)

# Load world cities data and keep only 300 cities, which will give us 90,000 pairs
data(world.cities)
world.cities <- world.cities[which(world.cities$country.etc=="Israel")[1:5],] 

# create all possible pairs of origin-destination in a long format
dt <- expand.grid.df(world.cities,world.cities)
names(dt)[10:11] <- c("lat_dest","long_dest")

# calculate distances in meters:
setDT(dt)[ , dist_km := distGeo(matrix(c(long, lat), ncol = 2), 
                                matrix(c(long_dest, lat_dest), ncol = 2))/1000]

  
calcDistExpand <- function(timegroupDF) {#this method gets the df section within a time group and exapands distance grids 
  library(reshape) ;library(data.table)
  #columns relevant: lat, ID, long
  
  }
  #

