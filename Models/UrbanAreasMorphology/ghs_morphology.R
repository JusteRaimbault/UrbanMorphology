
setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Models/UrbanAreasMorphology'))

library(rgdal)
#library(sf)
library(raster)
library(ggplot2)
library(reshape2)

source('../MorphologyPackage/morphology.R')

ucdb <- readOGR(paste0(Sys.getenv('CS_HOME'),'/Data/JRC_EC/GHS/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0'),'GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0')
#ucdb <- st_read(paste0(Sys.getenv('CS_HOME'),'/Data/JRC_EC/GHS/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0'),'GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0')

years = c('1975','1990','2000','2015')

pops = list()

for(year in years){
  pops[[year]] = raster(paste0(Sys.getenv('CS_HOME'),'/Data/JRC_EC/GHS/GHS_POP_GPW4',year,'_GLOBE_R2015A_54009_1k_v1_0/GHS_POP_GPW4',year,'_GLOBE_R2015A_54009_1k_v1_0.tif'))
}

largest_ua <- ucdb[which(ucdb@data[,'P15']==max(ucdb@data[,'P15'])),]
largest_ua <- spTransform(largest_ua,crs(pops[['2015']]))

currentua = spTransform(ucdb[order(ucdb@data[,'P15'],decreasing = T)[100],],crs(pops[['2015']])) #by chance Barcelona !
# currentua=largest_ua

currentpop = pops[['2015']]

#colFromX(pops[['2015']],100.0)
#cellsFromExtent(pops[['2015']],newExtent())
cells <- cellFromPolygon(currentpop,currentua)
rows <- rowFromCell(currentpop,cells[[1]])
cols <- colFromCell(currentpop,cells[[1]])

#xx <- xFromCell(currentpop,cells[[1]])
#yy <- yFromCell(currentpop,cells[[1]])

popvalues <- getValuesBlock(currentpop,row=min(rows)-(max(rows)-min(rows))/2,nrows=2*(max(rows)-min(rows)),col=min(cols)-(max(cols)-min(cols))/2,ncols=2*(max(cols)-min(cols)),format='matrix')

#g=ggplot(melt(popvalues),aes(x=Var1,y=Var2,fill=value))
#g+geom_raster()

#g=ggplot(melt(popvalues),aes(x=Var1,y=Var2,fill=value>100))
#g+geom_raster()

t = proc.time()
moranIndex(popvalues)
show(proc.time()[3] - t[3])

t = proc.time()
convolMoran(raster(popvalues))
show(proc.time()[3] - t[3])
# -> fastest
t = proc.time()
averageDistance(popvalues)
show(proc.time()[3] - t[3])


t = proc.time()
moran(popvalues)
show(proc.time()[3] - t[3])

# indicators ? ! come back to Le Nechet, 2015: add fractal dimension, acentrism.
# + pb of overlap if on too much areas ? take relative size only

acentrism(popvalues)
acentrism(popvalues,quantiles=seq(0.0,0.9,0.1))
acentrism(popvalues,quantiles=seq(0.0,0.95,0.05))
acentrism(popvalues,quantiles=seq(0.0,0.99,0.01))
acentrism(popvalues,quantiles=seq(0.0,0.999,0.001))
# somehow converges for the integral estimator - still performances issues ?

fractaldimension(popvalues)


