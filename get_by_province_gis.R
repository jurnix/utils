setwd("~/projects/rachel/thai gis/all")

library(rgdal) 
library(raster)
library(rgeos)


#Read in map file, adjust file name
#map=readOGR("THA_adm1.shp","THA_adm1")
#lake=map[map$ID_1==c(39,64),]
#map=map[map$ID_1[c(prov$mapid)],]
#map=spTransform(map, CRS("+proj=longlat"))
#coords=coordinates(map)

watercov=readOGR("THA_water_areas_dcw.shp","THA_water_areas_dcw")
water=readOGR("THA_water_lines_dcw.shp","THA_water_lines_dcw")
rail=readOGR("THA_rails.shp","THA_rails")
road=readOGR("THA_roads.shp","THA_roads")
cov=raster("THA_cov.grd")
alt=raster("THA_alt.grd")

provinces=readOGR("THA_adm1.shp", "THA_adm1")

#rowsName<-c('ID','Name','roads','rails','waterLines','waterCover','altitudeMean','altitudeStd', 'Centroid')
output<-data.frame()#rows.name=rowsName)

lids<-list()
lnames<-list()
lroads<-list()
lrails<-list()
larea<-list()
lwaterLines<-list()
lwaterCover<-list()
laltitudeMean<-list()
laltitudeStd<-list()
lCentroid<-list()
lcenterlong<-list()
lcenterlat<-list()

#########################################################################

for (i in 1:length(provinces) ) {
  #bar.squared[i] = bar[i]^2
  
  province<-subset(provinces, provinces$ID_1 == i)
  
  str(province$NAME_1)
  str(province$ID_1)
  
  altvalues<-extract(alt, province)
  altcrop<-crop(alt, province)
  
  # roads
  # get all roads found inside the province amnat
  roadsamnat<-gIntersection(road, province)
  # distancia total de línies en km
  if (!is.null(roadsamnat)){
    lengthroad<-SpatialLinesLengths(roadsamnat, longlat=TRUE)
  }else {
    lengthroad<-0.0
  }
  
  # rail lines
  railamnat<-gIntersection(rail, province)
  # check if there rail lines in region
  if (!is.null(railamnat)){
    lengthrail<-SpatialLinesLengths(railamnat, longlat=TRUE)
  } else{
    lengthrail<-0.0
  }
  
  # water coverage    
  waterbyindex<-which(gIntersects(watercov, province,byid=TRUE))
  # there is water in region?
  if (length(waterbyindex) > 0){
    watercovamnat<-polygons(watercov)[waterbyindex]
    watercovArea<-sum(areaPolygon(watercovamnat)) / 1000000 # get altitude mean and st dev
  }
  else {
    watercovArea=0.0
  }
  
  # water lines
  wateramnat<-gIntersection(water, province)
  if (!is.null(wateramnat)){
    lengthwater<-SpatialLinesLengths(wateramnat, longlat=TRUE)
  }else{
    lengthwater<-0.0
  }
  
  # get province area in m square
  amnatArea<-areaPolygon(province)/1000000# get altitude mean and st dev
  
  # altitude
  alt.mean <- lapply(altvalues, FUN=mean, na.rm=TRUE)
  alt.sd<- lapply(altvalues, FUN=sd, na.rm=TRUE)
  
  # calc coords center of polygon
  trueCentroids<- gCentroid(province, byid=TRUE)
    
  lids<-append(lids, province$ID_1)
  lnames<-append(lnames, levels(province$NAME_1)[[i]])
  lroads<-append(lroads, lengthroad)
  lrails<-append(lrails, lengthrail)
  lwaterLines<-append(lwaterLines, lengthwater)
  lwaterCover<-append(lwaterCover, watercovArea)
  laltitudeMean<-append(laltitudeMean, alt.mean[[1]])
  laltitudeStd<-append(laltitudeStd, alt.sd[[1]])  
  lcenterlong<-append(lcenterlong, trueCentroids@coords[[1]])
  lcenterlat<-append(lCentroid, trueCentroids@coords[[2]])
  larea<-append(larea, amnatArea)
  
  jpeg(paste(i, "_layers.jpg"))
  plot(altcrop)
  plot(province, add=T,border="red")
  if (!is.null(railamnat)){
    plot(roadsamnat,add=T, border="green")
  }
  if (!is.null(railamnat)){
    plot(railamnat,add=T, border="gray")
  }
  if (length(waterbyindex) > 0){
    plot(watercovamnat,add=T, border="blue")
  }
  if (!is.null(wateramnat)){
    plot(wateramnat,add=T, border="blue")
  }
  plot(trueCentroids,add=T)
  dev.off() 
}

rowsName<-c('ID','Name', 'area','roads','rails','waterLines','waterCover','altitudeMean','altitudeStd', 'cLatitude', 'cLongitude')
df<-data.frame(unlist(lids), unlist(lnames), unlist(larea), unlist(lroads), unlist(lrails), unlist(lwaterLines)
               , unlist(lwaterCover), unlist(laltitudeMean), unlist(laltitudeStd), unlist(lcenterlat), unlist(lcenterlong))
colnames(df) <- rowsName
write.csv(df, file="thailand_data.csv")
