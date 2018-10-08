
setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology'))

library(dplyr)
library(sp)

source('Models/StaticCorrelations/mapFunctions.R')

# network and pop indicators
areasize=100;offset=50;factor=0.5
indics = as.tbl(loadIndicatorData(paste0(Sys.getenv('CS_HOME'),"/UrbanMorphology/Data/StaticCorrelations/res/europecoupled_areasize",areasize,"_offset",offset,"_factor",factor,"_temp.RData")))

# emissions
# -> aggregate on different gases: for rows emissions -> indics, no need for now (same grid each gas)
emissions <- as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_C_2012/v432_CO2_excl_short-cycle_org_C_2012.txt',sep=';'))

# filter on bbox
emissions=emissions[emissions$lon>=min(indics$lonmin)&emissions$lon<=max(indics$lonmin)&emissions$lat>=min(indics$latmin)&emissions$lat<=max(indics$latmin),]

# by hand "join" (aggreg)

#library(doParallel)
#cl <- makeCluster(4,outfile='log')
#registerDoParallel(cl)
#rows <- foreach(i=1:nrow(emissions)) %dopar% {
rows=c();mins=c()
#for(i in 1:nrow(emissions)){
for(i in 1:nrow(indics)){
  #dif=abs(emissions$lon-indics$lonmin[i])+abs(emissions$lat-indics$latmin[i])
  #dif=abs(emissions$lon[i]-indics$lonmin)+abs(emissions$lat[i]-indics$latmin)
  if(i%%1000==0){show(paste0(100*i/nrow(indics),'%'))}
  #dif=sqrt((emissions$lon[i]-indics$lonmin)^2+(emissions$lat[i]-indics$latmin)^2)
  #dif = spDistsN1(as.matrix(indics[,c("lonmin","latmin")]),c(emissions$lon[i],emissions$lat[i]),longlat = T)
  dif = spDistsN1(as.matrix(emissions[,c("lon","lat")]),c(indics$lonmin[i],indics$latmin[i]),longlat = T)
  #indmin = ifelse(min(dif) < 15 ,which(dif==min(dif))[1],NA) # check distrib of distances
  indmin = which(dif==min(dif))[1]
  mins=append(mins,min(dif))
  rows=append(rows,indmin)
}
#stopCluster(cl)

#hist(log(mins),breaks=500)
#hist(mins[mins<20],breaks=500)
# need to check visually
#library(ggplot2)
#library(viridis)
#g=ggplot(data.frame(emissions,mins)[mins<20,],aes(x=lon,y=lat,fill=mins))
#g+geom_raster()+scale_fill_viridis(option = "magma", direction = -1)
#ggsave(file=paste0('Results/Maps/Test/spatialaliasing.png'),width=30,height=20,units='cm')
#length(which(mins<20)) # -> actually the resolution of emissions is smaller than our indics : redo with rows for indics

# keep all rows, put na in consolidation by keepping also mins
#save(rows,mins,file='Models/DataPreparation/rowsemissions.RData')
save(rows,mins,file='Models/DataPreparation/rowsindics.RData')
#system('rm log')
#stopCluster(cl)


