
setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Models/Percolation'))

library(dplyr)
#library(viridis)

source('../StaticCorrelations/mapFunctions.R')
source('percolationFunctions.R')

areasize=100;offset=50;factor=0.5
indics = as.tbl(loadIndicatorData(paste0(Sys.getenv('CS_HOME'),"/UrbanMorphology/Data/StaticCorrelations/res/europecoupled_areasize",areasize,"_offset",offset,"_factor",factor,"_temp.RData")))

# TODO add emission loading (transportation ?) and spatial join

quantiles = c(0.75,0.8,0.85,0.9,0.95)
radiuses=c(8000,10000,15000,20000,50000,75000,100000)
nwindics= c("ecount","mu","vcount")

params=matrix(0,length(quantiles)*(length(quantiles)+1)*length(radiuses)*length(nwindics),4)
i=1
for(nwcol in nwindics){for(nwthq in c(0,quantiles)){for(popthq in quantiles){
  for(radius in radiuses){
    params[i,]=c(nwcol,popthq,nwthq,radius)
    i=i+1
  }
}}}
params=as.data.frame(params)
colnames(params)=c("nwcol","popthq","nwthq","radius")
params$nwthq=as.numeric(as.character(params$nwthq));params$popthq=as.numeric(as.character(params$popthq));params$radius=as.numeric(as.character(params$radius))

library(doParallel)
cl <- makeCluster(50,outfile='log')
registerDoParallel(cl)
res <- foreach(i=1:nrow(params)) %dopar% {
  source('percolationFunctions.R')
  #p = conditionalPercolation(d=indics,popthq=params$popthq[i],nwthq=params$nwthq[i],radius=params$radius[i])
  return(conditionalPercolation(d=indics,popthq=params$popthq[i],nwthq=params$nwthq[i],radius=params$radius[i]))
}
stopCluster(cl)

save(res,file='directSampling.RData')






