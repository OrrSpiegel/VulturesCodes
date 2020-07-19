## this function is called by the main one ReadMoveBankDataFormat
## it loops on all time groups
## for each one it updates the nubers of co-occuring vultures, including self 
## also the updates the counters of simulatnounsly tracked vultures 
## note that each dyad, including self is counter twice A--A and B--A or if self A--A and A--A
#takes ~3.5 hours

ColumToSelect=c("ID","location_lat","location_long","Easting","Northing","timegroup","group")
  
for (timgrpind in 1: max(DatasetOhadF$timegroup)){#loop on all time groups
  ## extract current time group (#18458 has a good example)
  #subset(DatasetOhadF, timegroup==timgrpind,select=c("ID","location_lat","location_long"))
  timegroupDF=subset(DatasetOhadF, timegroup==timgrpind,select=ColumToSelect)
  #print(timegroupDF[order(timegroupDF$group),])#,'timegroupDF$ID'#plot it
  
  ## working within this time group: dyads and distances: 
  timegroupDF=subset(timegroupDF, timegroup==timgrpind,select=c("ID","location_lat","location_long","group"))
  DT <- expand.grid.df(timegroupDF,timegroupDF);  names(DT)[5:7] <- c('ID2',"lat_secondInd","long_secondInd")
  #distances:
  setDT(DT)[ , dist_km := distGeo(matrix(c(location_long, location_lat), ncol = 2), 
                                  matrix(c(long_secondInd, lat_secondInd), ncol = 2))/1000];
  
  ## create all possible dyads of vultures in a long format:
  PresentVultures=subset(DT, dist_km==0&as.character(ID)==as.character(ID2),select=c(1,5))#these are the vultures present in this time group
  PresentVultures=as.data.frame(unique(PresentVultures$ID))#PresentVultures=unique(PresentVultures$ID);
  PresentVultures=expand.grid.df(PresentVultures,PresentVultures,unique=F)#these are the dyads that has concurrent time point
  names(PresentVultures)=c('ID','ID2')
  
  ## loop on current dyads (including self) to update co-occurances:
  for (dyadcnt in 1: dim(PresentVultures)[1]){#length just half since each dyad is counted twice AB and BA
    Dyadind=which(SimlDataPntCnt$ID==PresentVultures$ID[dyadcnt]&SimlDataPntCnt$ID2==PresentVultures$ID2[dyadcnt])
    SimlDataPntCnt$counter[Dyadind]=SimlDataPntCnt$counter[Dyadind]+1;#
  }#for loop on current dyads 
  
  #if(CountTwice){
  ## since self dyads appear only once:
  SelfDyad=which(PresentVultures$ID==PresentVultures$ID2);#since self dyads appear only once, another loop on them, so diag will be counted twice like the rest
  for (dyadcnt in (SelfDyad)){#count dyad twice? A-B and B-A
    Dyadind=which(SimlDataPntCnt$ID==PresentVultures$ID[dyadcnt]&SimlDataPntCnt$ID2==PresentVultures$ID2[dyadcnt])
    SimlDataPntCnt$counter[Dyadind]=SimlDataPntCnt$counter[Dyadind]+1;#
  }#count dyad twice? A-B and B-A
  #so the diagonal, self is counted twice since its A-A and A-A later. others are A-B and B-A, 
  
  ## now setting interacting dyads
  #if (dim(timegroupDF)[1]>=1){#more than one vulture in this time group?
  InteractingSelf=subset(DT, dist_km==0 & (as.character(ID)==as.character(ID2))) # just including self interactions, once, not multiple times
  InteractingSelf=InteractingSelf[!duplicated(InteractingSelf$ID),]# just including self interactions, once, not multiple times
  InteractingDyads=subset(DT,(dist_km<=DistThresholM/1000 & (as.character(ID)!=as.character(ID2)))) # not including self interactions
  InteractingDyads=InteractingDyads[!duplicated(InteractingDyads[,c("ID","ID2")])];
  
  if(dim(PresentVultures)[1]<dim(InteractingDyads)[1])      {
    break
    }#for debugging
  
  InteractingDyads=rbind(InteractingDyads,InteractingSelf,InteractingSelf)#counting the self twice to keep it like the twice of the non-self dyads
  ## a loop on interacting dyads, non-self, for updating the CoocurCountr storge
  for (dyadcnt in 1: dim(InteractingDyads)[1]){#
    Dyadind=which(CoocurCountr$ID==InteractingDyads$ID[dyadcnt]&CoocurCountr$ID2==InteractingDyads$ID2[dyadcnt])
    CoocurCountr$counter[Dyadind]=CoocurCountr$counter[Dyadind]+1;#
  }#loop 
  
  # }#if making sure Ids match
  
  ## debugging:
  #print(InteractingDyads[order(InteractingDyads$group,InteractingDyads$ID),])#,'timegroupDF$ID'#plot it
  #conditions for debugging:
  #(dim(InteractingSelf)[1]*2<dim(InteractingDyads)[1])
  #if(length(unique(PresentVultures$ID))>1)
  #if(dim(InteractingSelf)[1]*2<dim(InteractingDyads)[1]){#read dyads interact
  #  print(timgrpind)
  #  print('timgrpind')
  #}
  
  rm(list=c("DT","InteractingDyads","InteractingSelf","SelfDyad","Dyadind","dyadcnt","timegroupDF"))
  
}#loop on time groups

