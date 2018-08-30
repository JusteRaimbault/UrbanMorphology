
setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology'))

library(dplyr)
library(ggplot2)
library(raster)
library(viridis)
library(rgdal)
library(rgeos)

source('Models/StaticCorrelations/mapFunctions.R')


#indics <-loadIndicatorData('Data/StaticCorrelations/res/20150806_europe50km_10kmoffset_100x100grid.csv')
areasize=100;offset=50;factor=0.5
indics = as.tbl(loadIndicatorData(paste0("Data/StaticCorrelations/res/europecoupled_areasize",areasize,"_offset",offset,"_factor",factor,"_temp.RData")))
#indics=indics[,1:10]
indics=indics[apply(indics,1,function(r){prod(as.numeric(!is.na(r)))>0}),]

emissions <- as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_C_2012/v432_CO2_excl_short-cycle_org_C_2012.txt',sep=';'))

# load spatial mask to select area
countrycode='FR'
countries = readOGR('Data/gis','countries')
country = countries[countries$CNTR_ID==countrycode,]
datapoints = SpatialPoints(data.frame(emissions[,c("lon","lat")]),proj4string = countries@proj4string)
selectedpoints = gContains(country,datapoints,byid = TRUE)
sdata = emissions[selectedpoints,]

datapoints = SpatialPoints(data.frame(indics[,c("lonmin","latmin")]),proj4string = countries@proj4string)
selectedpoints = gContains(country,datapoints,byid = TRUE)
sindics = indics[selectedpoints,]

#g=ggplot(indics,aes(x=lonmin,y=latmin,fill=moran))
#g+geom_raster()+scale_fill_viridis(option = "magma", direction = -1)
#ggsave(file=paste0('Results/Maps/Test/moran.png'),width=30,height=20,units='cm')

#g=ggplot(sdata,aes(x=lon,y=lat,fill=norm_log_emission))
#g+geom_raster()+scale_fill_viridis(option = "magma", direction = -1)
#ggsave(file=paste0('Results/Maps/Test/CO2_2012_FR.png'),width=30,height=20,units='cm')


#library(doParallel)
#cl <- makeCluster(60,outfile='log')
#registerDoParallel(cl)
#rows <- foreach(i=1:nrow(indics)) %dopar% {
#  dif=abs(emissions$lon-indics$lonmin[i])+abs(emissions$lat-indics$latmin[i])
#  return(which(dif==min(dif)))
#}
#stopCluster(cl)

smootheddata=data.frame()
for(i in 1:nrow(sindics)){
  dif=abs(sdata$lon-sindics$lonmin[i])+abs(sdata$lat-sindics$latmin[i])
  smootheddata=rbind(smootheddata,colMeans(sdata[which(dif<quantile(dif,0.05)),]))
}

# -> only correlated to pop in highly populated areas ? -> smoothing may increase significantly
cor(log(sdata[rows,c(-1,-2)]),log(sindics[,c(-1,-2)]))

cor(log(smootheddata[,c(-1,-2)]),log(sindics[,c(-1,-2)]))
cor(smootheddata[,c(-1,-2)],sindics[,c(-1,-2)])

library(corrplot)
corrplot(cor(smootheddata[,c(-1,-2)],sindics[,c(-1,-2)]),method="color")

#g=ggplot(sdata[rows,],aes(x=lon,y=lat,fill=norm_log_emission))
#g+geom_raster()+scale_fill_viridis(option = "magma", direction = -1)
#g=ggplot(sindics,aes(x=lonmin,y=latmin,fill=log(totalPop)))
#g+geom_raster()+scale_fill_viridis(option = "magma", direction = -1)







