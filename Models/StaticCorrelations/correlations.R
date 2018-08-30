
setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology'))

library(dplyr)
library(ggplot2)
library(raster)
library(viridis)

source('Models/StaticCorrelations/mapFunctions.R')


#indics <-loadIndicatorData('Data/StaticCorrelations/res/20150806_europe50km_10kmoffset_100x100grid.csv')
areasize=100;offset=50;factor=0.5
indics = as.tbl(loadIndicatorData(paste0("Data/StaticCorrelations/res/europecoupled_areasize",areasize,"_offset",offset,"_factor",factor,"_temp.RData")))
indics=indics[,1:10]
indics=indics[apply(indics,1,function(r){prod(as.numeric(!is.na(r)))>0}),]

emissions <- as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_C_2012/v432_CO2_excl_short-cycle_org_C_2012.txt',sep=';'))

#g=ggplot(indics,aes(x=lonmin,y=latmin,fill=moran))
#g+geom_raster()+scale_fill_viridis(option = "magma", direction = -1)
#ggsave(file=paste0('Results/Maps/Test/moran.png'),width=30,height=20,units='cm')

library(doParallel)
cl <- makeCluster(60,outfile='log')
registerDoParallel(cl)

rows <- foreach(i=1:nrow(indics)) %dopar% {
  dif=abs(emissions$lon-indics$lonmin[i])+abs(emissions$lat-indics$latmin[i])
  return(which(dif==min(dif)))
}
save(rows,file='Models/StaticCorrelations/rowstmp.RData')
system('rm log')

stopCluster(cl)
