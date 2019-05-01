# Apr 23rd :-) 2019 trying to read movebank data using move

require(move)
#require(dplyr) 
#require(tidyverse) 
require(ggmap) #these packages are necessary to work with google maps
require(mapproj)

#movebank user ors PW: eesH1eeF
MB.LoginObject=movebankLogin(username='ors',password='eesH1eeF')
#"Gyps fulvus INPA Hatzofe" Movebank ID	6071688
#"HUJ MoveEcol Lab Israel: Griffon vulture Gyps fulvus" Movebank ID	6638215
#"E-obs Gsm Vultures Israel" Movebank ID	7359070
#timestamp_start, timestamp_end character or POSIXct=’yyyyMMddHHmmssSSS’
MoveStackDatasetOhad=getMovebankData(study=6071688, login=MB.LoginObject,
        includeExtraSensors=FALSE, deploymentAsIndividuals=FALSE,removeDuplicatedTimestamps=TRUE)#animalName
MoveStackDatasetOrr=getMovebankData(study=6638215, login=MB.LoginObject,
                                 includeExtraSensors=FALSE, deploymentAsIndividuals=FALSE,removeDuplicatedTimestamps=TRUE)#animalName

#returns or a MoveStack object can use as.data.frame
#corridor
show(MoveStackDatasetOhad)
summary(MoveStackDatasetOrr)
citations(MoveStackDatasetOrr)
equalProj(MoveStackDatasetOrr)
n.locs(MoveStackDatasetOhad)
timeLag(MoveStackDatasetOhad)


DatasetOrr=as.data.frame(MoveStackDatasetOrr) #converting to a dataset 
DatasetOhad=as.data.frame(MoveStackDatasetOhad) #converting to a dataset 
#str(DatasetOrr) ;str(DatasetOhad)
# removing un-necessary columns
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
                                          'acceleration_raw_x','acceleration_raw_z',"acceleration_raw_z","eobs_horizontal_accuracy_estimate","eobs_speed_accuracy_estimate");  DatasetOhad=DatasetOhad[!VarsToRemove]


#DatasetOrr$update_ts=NULL;DatasetOrr$taxon_canonical_name=DatasetOrr$TBA=NULL #Orr

#myvars <- names(mydata) %in% c("v1", "v2", "v3")
#unique(DatasetOrr$sensor_type_id)
unique(DatasetOhad$gps_satellite_count)
View(DatasetOhad)
DatasetOhadF$height_above_ellipsoid
### nitial dataset exploration & simple filteres ######
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


#another try without going to a data frame since i cannot go back
unstackedOhad=split(MoveStackDatasetOhad)
#summary(unstackedOhad)  
PercentThrown=vector(length=length(unstackedOhad))
for (indv in 1:length(unstackedOhad) ){#loop on individuals, now in separate Move objects
  plot(unstackedOhad[[indv]])
  dim(unstackedOhad[[indv]]@data)
  dim(unstackedOhad[[indv]]@coords)
  
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
  dim(unstackedOhad[[indv]]@data)
  #choosing indices to keep for this individual
  
  indx=1:dim(unstackedOhad[[indv]]@data)[1]
  if(sum(unstackedOhad[[indv]]@data$heading <= 360,na.rm=T)){
    indx=intersect(indx,which(unstackedOhad[[indv]]@data$heading <= 360))}
  if(sum(unstackedOhad[[indv]]@data$ground_speed<=120,na.rm=T)){
    indx=intersect(indx,which(unstackedOhad[[indv]]@data$ground_speed<=120))}
  #if(sum(unstackedOhad[[indv]]@data$gps_time_to_fix<=89,na.rm=T)){
  #   indx=intersect(indx,which(unstackedOhad[[indv]]@data$gps_time_to_fix<=89))}
  if(sum(unstackedOhad[[indv]]@data$gps_satellite_count>=3,na.rm=T)){
    indx=intersect(indx,which(unstackedOhad[[indv]]$gps_satellite_count>=3))}
  
  ###subsetting the different slots of this move object
  print(paste('i throw out', dim(unstackedOhad[[indv]]@data)[1]-length(indx), 'points, out of',dim(unstackedOhad[[indv]]@data)[1]))
  PercentThrown[indv]=(dim(unstackedOhad[[indv]]@data)[1]-length(indx))/dim(unstackedOhad[[indv]]@data)[1]
  unstackedOhad[[indv]]@data=unstackedOhad[[indv]]@data[indx,]
  unstackedOhad[[indv]]@coords=unstackedOhad[[indv]]@coords[indx,]
  #unstackedOhad[[indv]]@idData
  #unstackedOhad[[indv]]@proj4string
  unstackedOhad[[indv]]@sensor=unstackedOhad[[indv]]@sensor[indx]
  unstackedOhad[[indv]]@timestamps=unstackedOhad[[indv]]@timestamps[indx]
  
  lines(unstackedOhad[[indv]],col='red')
  #plot(unstackedOhad[[indv]], type="o", col=3, lwd=2, pch=20, xlab="location_long", ylab="location_lat")
  
  #plotting a google map
  DF1indiv=as.data.frame(unstackedOhad[[indv]])
  m <- get_map(bbox(extent(min(DF1indiv$location_long), 
                           max(DF1indiv$location_long), 
                           min(DF1indiv$location_lat), 
                           max(DF1indiv$location_lat) )*1.1), source="stamen", zoom=10)
  ggmap(m)+geom_path(data=DF1indiv, aes(x=location_long, y=location_lat))

  
  }#loop


#plotting a google map ####
for (indv in 1:length(unstackedOhad) ){#loop on individuals, now in separate Move objects
  DF1indiv=as.data.frame(unstackedOhad[[indv]])
  m <- get_map(bbox(extent(min(DF1indiv$location_long), 
                           max(DF1indiv$location_long), 
                           min(DF1indiv$location_lat), 
                           max(DF1indiv$location_lat) )*1.1), source="stamen") #can be zoom 8 for guys with LFRs and 12 or local guys
  ggmap(m)+geom_path(data=DF1indiv, aes(x=location_long, y=location_lat))
}#loop on individuals




#### Calculation of the utilization distribution ####
brownian.bridge.dyn