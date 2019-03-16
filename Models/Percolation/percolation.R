
setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Models/Percolation'))

library(dplyr)

#source('Models/Percolation/percolationFunctions.R')
source('percolationFunctions.R')

# assumes data has been consolidated before
load('../../Data/consolidated/indics.RData')

#purpose='directSampling'
purpose='test'

# parameter values
#popquantiles=c(0.85,0.9,0.95)
popquantiles=c(0.95)

#nwquantiles=c(0.0,0.8,0.9,0.95)
nwquantiles=c(0.95)

#radiuses=c(8000,10000,15000,20000,50000)
radiuses=c(8000,10000)

#nwindics= c("ecount","mu","vcount","euclPerf")
nwindics= c("ecount")

#gammas=c(0.5,1,1.5,2.0)
gammas=c(1)
#decays=c(100,1000,10000,50000,100000)
decays=c(1000)

params=matrix(0,length(popquantiles)*length(nwquantiles)*length(radiuses)*length(nwindics)*length(gammas)*length(decays),6)
i=1
for(nwcol in nwindics){for(nwthq in nwquantiles){for(popthq in popquantiles){for(radius in radiuses){for(gamma in gammas){for(decay in decays){
    params[i,]=c(nwcol,popthq,nwthq,radius,gamma,decay)
    i=i+1
  }}}}}}
params=as.data.frame(params)
colnames(params)=c("nwcol","popthq","nwthq","radius","gamma","decay")
params$nwthq=as.numeric(as.character(params$nwthq));params$popthq=as.numeric(as.character(params$popthq));params$radius=as.numeric(as.character(params$radius));params$gamma=as.numeric(as.character(params$gamma));params$decay=as.numeric(as.character(params$decay))

show(paste0('Number of parameters points = ',nrow(params)))

library(doParallel)
#cl <- makeCluster(60,outfile='log')
cl <- makeCluster(2,outfile='log')

registerDoParallel(cl)
res <- foreach(i=1:nrow(params)) %dopar% {
#res=list()
#for(i in 1:nrow(params)){
    #source('Models/Percolation/percolationFunctions.R')
  source('percolationFunctions.R')
  # Mandatory args : d,radius,popthq,nwcol,nwthq,gamma,decay
  return(graphPercolation(d=indics,
                                radius=params$radius[i],
                                popthq=params$popthq[i],
                                nwcol=as.character(params$nwcol[i]),
                                nwthq=params$nwthq[i],
                                gamma=params$gamma[i],
                                decay=params$decay[i]
                                )
    )
}
stopCluster(cl)

save(res,file=paste0('res/',purpose,'_',format(Sys.time(), "%Y%m%d_%H%M%S"),'.RData'))

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







