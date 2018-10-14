
setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Models/Percolation'))

library(dplyr)

#source('Models/Percolation/percolationFunctions.R')
source('percolationFunctions.R')

# assumes data has been consolidated before
load('../../Data/consolidated/indics.RData')

# parameter values
quantiles = c(0.75,0.8,0.85,0.9,0.95)
radiuses=c(8000,10000,15000,20000,50000,75000,100000)
nwindics= c("ecount","mu","vcount","euclPerf")
gammas=c(0.5,1,1.5)
decays=c(10,100,500,1000)

params=matrix(0,length(quantiles)*(length(quantiles)+1)*length(radiuses)*length(nwindics)*length(gammas)*length(decays),6)
i=1
for(nwcol in nwindics){for(nwthq in c(0,quantiles)){for(popthq in quantiles){for(radius in radiuses){for(gamma in gammas){for(decay in decays){
    params[i,]=c(nwcol,popthq,nwthq,radius,gamma,decay)
    i=i+1
  }}}}}}
params=as.data.frame(params)
colnames(params)=c("nwcol","popthq","nwthq","radius","gamma","decay")
params$nwthq=as.numeric(as.character(params$nwthq));params$popthq=as.numeric(as.character(params$popthq));params$radius=as.numeric(as.character(params$radius));params$gamma=as.numeric(as.character(params$gamma));params$decay=as.numeric(as.character(params$decay))

library(doParallel)
cl <- makeCluster(60,outfile='log')
registerDoParallel(cl)
res <- foreach(i=1:nrow(params)) %dopar% {
  #source('Models/Percolation/percolationFunctions.R')
  source('percolationFunctions.R')
  #p = conditionalPercolation(d=indics,popthq=params$popthq[i],nwthq=params$nwthq[i],radius=params$radius[i])
  # Mandatory args : d,radius,popthq,nwcol,nwthq,gamma,decay
  return(aggregIndics(conditionalPercolation(d=indics,
                                radius=params$radius[i],
                                popthq=params$popthq[i],
                                nwcol=as.character(params$nwcol[i]),
                                nwthq=params$nwthq[i],
                                gamma=params$gamma[i],
                                decay=params$decay[i]
                                )
         ))
}
stopCluster(cl)

save(res,file=paste0('res/directSampling_',format(Sys.time(), "%Y%m%d_%H%M%S"),'.RData'))

#############

# load('directSampling.RData')

# get summary indicators of pop / areas distrib

# areasalpha = unlist(lapply(res,function(l){
#   areas = l$areas
#   return(ifelse(length(which(areas>0)),lm(data=data.frame(rank=log(1:length(areas[areas>0])),size=log(areas[areas>0])),size~rank)$coefficients[2],NA))
# }))
# 
# popsalpha = unlist(lapply(res,function(l){
#   pops = l$pops
#   return(ifelse(length(which(pops>0)),lm(data=data.frame(rank=log(1:length(pops[pops>0])),size=log(pops[pops>0])),size~rank)$coefficients[2],NA))
# }))
# 
# clustnum = unlist(lapply(res,function(l){length(l$pops)}))
# 
# d=data.frame(areasalpha,popsalpha,clustnum)
# 
# g=ggplot(d[clustnum>15,],aes(x=areasalpha,y=popsalpha,colour=clustnum))
# g+geom_point()
# 
# # -> hierarchy of "endogenous city system" ?
# # : find criteria on which to optimize !
# 







