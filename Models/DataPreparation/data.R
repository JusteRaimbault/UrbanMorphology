
setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology'))

library(dplyr)

source('Models/StaticCorrelations/mapFunctions.R')

# network and pop indicators
areasize=100;offset=50;factor=0.5
indics = as.tbl(loadIndicatorData(paste0(Sys.getenv('CS_HOME'),"/UrbanMorphology/Data/StaticCorrelations/res/europecoupled_areasize",areasize,"_offset",offset,"_factor",factor,"_temp.RData")))

# emissions
# -> aggregate on different gases: for rows emissions -> indics, no need for now (same grid each gas)
emissions <- as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_C_2012/v432_CO2_excl_short-cycle_org_C_2012.txt',sep=';'))


# by hand "join" (aggreg)

library(doParallel)
cl <- makeCluster(10,outfile='log')
registerDoParallel(cl)
rows <- foreach(i=1:nrow(emissions)) %dopar% {
  #dif=abs(emissions$lon-indics$lonmin[i])+abs(emissions$lat-indics$latmin[i])
  #dif=abs(emissions$lon[i]-indics$lonmin)+abs(emissions$lat[i]-indics$latmin)
  show(paste0(100*i/nrow(emissions),'%'))
  dif=sqrt((emissions$lon[i]-indics$lonmin)^2+(emissions$lat[i]-indics$latmin)^2)
  indmin = which(dif==min(dif))[1]
  show(indmin)
  return(indmin)
}
stopCluster(cl)
save(rows,file='Models/DataPreparation/rowsemissions.RData')
system('rm log')
stopCluster(cl)


