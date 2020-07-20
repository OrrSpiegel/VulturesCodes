#this script  reads the kml file of roosts, finds the centroids and then claculate distances between them 
# clean workspace
rm(list=ls()) 
graphics.off()

## important parameters #####
 EdgeThrshldDist=40; #in kms

## load packages ##### 
require(rgdal)
require(igraph)

## read kml of roosts #####
#FileName='RoostsListAug2010.kml'
#LayerName='Roosting'
#Roosts = readOGR(dsn=FileName, layer=)

load('C:/Users/Noa/Dropbox/grants/NSF-BSF Orr/Net_fig/DataRoostForNoa1.RData')


# matrix of all the lsites and distances between them
RoostDistMtrxKm

## upload bipirtite network for oth social and interlayer links...





class(Roosts)#SpatialPolygonsDataFrame
summary(Roosts)#SpatialPolygonsDataFrame


## extract centroids and calculate the distance matrix, afer removing  one value #####
#RoostsNames=Roosts@data$Name
centroids <- getSpPPolygonsLabptSlots(Roosts)
RoostsCentroids=data.frame(centroids);RoostsCentroids$Names=Roosts@data$Name;names(RoostsCentroids)=c('Long', 'Lat','Name')
RoostsCentroids=subset(RoostsCentroids, Name != 'Other',select=c(Name,  Long, Lat)) 


#removing three far locations in saudi Arabia
RoostsCentroids=subset(RoostsCentroids, Name != 'SA_TabukNE') #
RoostsCentroids=subset(RoostsCentroids, Name != 'SA_Mughayra')
RoostsCentroids=subset(RoostsCentroids, Name != 'RamUmEsreenJordan')

## plot 
plot(Roosts, col = 'red', pbg="white")
points(centroids, pch = 3, col = "black")

## distance matrix in km between roosts
#RoostDistMtrxKm=dist(RoostsCentroids, diag=T, upper=T)
RoostDistMtrxKm <- raster::pointDistance(p1=RoostsCentroids[,c('Long','Lat')], lonlat=TRUE,allpairs=TRUE)/1000; #divided for kms!
RoostDistMtrxKm[upper.tri(RoostDistMtrxKm)] <- t(RoostDistMtrxKm)[upper.tri(RoostDistMtrxKm)]#filling the upper triangle
rownames(RoostDistMtrxKm)=RoostsCentroids$Name
colnames(RoostDistMtrxKm)=RoostsCentroids$Name


## checking within threshold distance otherwise exclude #####
RoostDistMtrxKm[RoostDistMtrxKm>EdgeThrshldDist] <- NA# if with NA then all lines are connected

## first plot of the matrix
image(RoostDistMtrxKm)
image(RoostDistMtrxKm, col=topo.colors(10),axes=T)#of heat colors or rainbow
title(main="RoostDistMtrxKm for roosts  \n ", font.main=4)
#colorRamp

## working with the Roosts Network ####
RoostDistMtrxKm[is.na(RoostDistMtrxKm)] <- 0# if with NA then all lines are connected
Graph2 <- graph.adjacency(adjmatrix=RoostDistMtrxKm, mode='undirected', diag = FALSE, weighted = TRUE)
plot.igraph(Graph2)
tkplot(Graph2,vertex.color = "yellow", vertex.frame.color = "gray20")

#creating a standardized location for vertexes
RoostsCentroids$GrphPos1=RoostsCentroids$Long-min(RoostsCentroids$Long);RoostsCentroids$GrphPos1=2*(RoostsCentroids$GrphPos1/max(RoostsCentroids$GrphPos1)-0.5)
RoostsCentroids$GrphPos2=RoostsCentroids$Lat -min(RoostsCentroids$Lat) ;RoostsCentroids$GrphPos2=2*(RoostsCentroids$GrphPos2/max(RoostsCentroids$GrphPos2)-0.5)
#range(RoostsCentroids$GrphPos1)
#plot(RoostsCentroids$GrphPos1,RoostsCentroids$Long)
#plot(RoostsCentroids$GrphPos2,RoostsCentroids$Lat)


plot.igraph(Graph2, #vertex.size = ideg*25 + 40, vertex.size2 = 30,
            vertex.color = "yellow", vertex.frame.color = "gray20",
            #vertex.shape = "rectangle", 
            #edge.arrow.size=0.5, edge.color='blue',
            layout = as.matrix(RoostsCentroids[,c('GrphPos1','GrphPos2')])
            )

tkplot(Graph2, #vertex.size = ideg*25 + 40, vertex.size2 = 30,
       vertex.color = "yellow", vertex.frame.color = "gray20",
       #vertex.shape = "rectangle", 
       #edge.arrow.size=0.5, edge.color='blue',
       layout = as.matrix(RoostsCentroids[,c('GrphPos1','GrphPos2')]))

########################################
######  ADDED BY NOA FOR NET PLOTTING ##################
########################################
ind_loc_for_bip=read.csv('C:/Users/Noa/Dropbox/grants/NSF-BSF Orr/Net_fig/ind_loc_bipartite.csv')


library('bipartite') # loaded here and not above cause otherwise it maskes some igraph function cause this package is based on sna...


### plot for grant proposal:###################
windows()
plot(Graph2,  vertex.color = 'tan3', vertex.size = 5,#vertex.label=NA,
		vertex.frame.color = "tan4", main='Distance between roosting sites')#, layout = as.matrix(RoostsCentroids[,c('GrphPos1','GrphPos2')]))

cfg1 <- fastgreedy.community(Graph2)  # dedect comunities

# plot the network with cluster colors:
# to make the colors match the figure in the paper:
cols1=gsub(4, 'tan4',cfg1$membership)
cols1=gsub(2, 'burlywood4',cols1)
cols1=gsub(3, 'chocolate3',cols1)
cols1=gsub(1, 'darkorange4',cols1)
cols1=gsub(5, 'darkorange3',cols1)

windows()
plot(Graph2,vertex.color=cols1, vertex.label=NA, mark.groups=communities(cfg1), vertex.size = 7, mark.col=c('tan2','burlywood2', 'tan','sandybrown', 'peachpuff'), 
		mark.border="orangered4", main='Distance between roosting sites',edge.color='orangered4')


windows()
plotweb(t(ind_loc_for_bip[,2:dim(ind_loc_for_bip)[2]]>5), method='normal',col.high='yellowgreen',bor.col.high='yellowgreen', 
		high.lablength=0,low.lablength=0,				
		col.low='tan4',bor.col.low='tan4',col.interaction='darkorange2',bor.col.interaction='darkorange2')

### create interactions networks
library('asnipe')
for_GBI=t(ind_loc_for_bip[,2:dim(ind_loc_for_bip)[2]])
for_GBI=as.data.frame(for_GBI)
names(for_GBI)=ind_loc_for_bip$X
net_gr=get_network(for_GBI>5, 
		data_format = "GBI",association_index = "SRI")

windows()
net1=graph.adjacency(as.matrix(net_gr),mode="undirected", weighted = TRUE, diag=FALSE)
plot(net1, vertex.label=NA, vertex.size = 8,vertex.color = 'yellowgreen', main='Roosting interactions')


cfg <- fastgreedy.community(net1)
# or: cfg <-cluster_fast_greedy(net)

# look at results:
cfg$membership

# plot the network with cluster colors:
# to make the colors match the figure in the paper:
cols=gsub(1, 'yellowgreen',cfg$membership)
cols=gsub(2, 'olivedrab',cols)


windows()
plot(net1,vertex.color=cols, vertex.label=NA, mark.groups=communities(cfg),mark.col=c('darkolivegreen1','darkolivegreen3'), 
		mark.border="black",edge.color='darkgreen', main='Roosting interactions')



rossting_net_mat=get.adjacency(net1)

rosts_net_mat=get.adjacency(Graph2)
save('rosts_net_mat','rossting_net_mat', 'net1', 'Graph2', file='C:/Users/Noa/Dropbox/grants/NSF-BSF Orr/Net_fig/nets.RData')



##############################################
####################### DRAFTS ############################################

cents <- coordinates(Roosts)
cents <- SpatialPointsDataFrame(coords=cents, data=Roosts@data,   proj4string=CRS("+proj=utm +zone=10+datum=WGS84"))
points(cents, col = "Blue")
writeSpatialShape(cents, "cents")



#geojsonio::centroid(Roosts)
centroids <- getSpPPolygonsLabptSlots(Roosts)
getPolygon(Roosts)
#Projection:
  proj4string(Roosts)
 proj4string(x) <-CRS("+proj=utm +zone=10+datum=WGS84")
 utmN <- '+proj=utm +zone=36        +ellps=WGS84 +datum=WGS84 +units=m +no_defs'  #north
 
 ## converting to UTM36North, but note that not all points are within, just the majority
 DatasetOhadF_utmN <- spTransform(DatasetOhadF_wgs, CRS(utmN))
 