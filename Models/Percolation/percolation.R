
setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Models/Percolation'))

library(dplyr)
#library(viridis)

source('../StaticCorrelations/mapFunctions.R')
source('percolationFunctions.R')

areasize=100;offset=50;factor=0.5
indics = as.tbl(loadIndicatorData(paste0(Sys.getenv('CS_HOME'),"/UrbanMorphology/Data/StaticCorrelations/res/europecoupled_areasize",areasize,"_offset",offset,"_factor",factor,"_temp.RData")))

quantiles = c(0.8,0.9,0.95,0.96,0.97,0.98,0.99)
radiuses=c(1000,2000,5000,10000,20000,50000,100000)
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
colnames(params)=c("nwcol","nwthq","popthq","radius")
params$nwthq=as.numeric(as.character(params$nwthq));params$popthq=as.numeric(as.character(params$popthq));params$radius=as.numeric(as.character(params$radius))

library(doParallel)
cl <- makeCluster(50,outfile='log')
registerDoParallel(cl)
res <- foreach(i=1:nrow(params)) %dopar% {
  source('percolationFunctions.R')
  return(conditionalPercolation(indics,params$popthq[i],params$nwthq[i],params$radius[i]))
}
stopCluster(cl)

save(res,file='directSampling.RData')






