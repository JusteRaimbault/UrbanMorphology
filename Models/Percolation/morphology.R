
# morphological indices

# TODO creating R packages
# https://cran.r-project.org/doc/contrib/Leisch-CreatingPackages.pdf

#setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Models/MorphologyPackage'))

library(units)

#######

library(raster)
library(sf)
library(Matrix)

#########
## Vector functions

#'
#'
moranPoints <-function(points,varcol,dmat=NULL){
  show(paste0('Moran for ',nrow(points),' points'))
  if(nrow(points)<5000){sppoints = points}else{
    th = quantile(points[[varcol]],c(1 - 5000/nrow(points)))
    sppoints = points[points[[varcol]]>th,]
  }
  if(is.null(dmat)){weights=spatialWeightsPoints(sppoints)}else{weights = Matrix(1/dmat);diag(weights)<-0}
  v = sppoints[[varcol]]
  v = v-mean(v)
  num = sum(Diagonal(x = v)%*%weights%*%Diagonal(x = v))
  denom = sum(v^2)
  normalization = nrow(points) / sum(weights)
  return(normalization*num/denom)
}

#'
#'
spatialWeightsPoints<-function(points){
  distances <- 1/st_distance(points)
  #distances[as.matrix(distances)==matrix(rep(Inf,nrow(distances)*ncol(distances)),nrow=nrow(distances))]=0
  # no point at the same place
  diag(distances)<-0
  mat = drop_units(distances)
  #mat[mat<0.01]=0
  #return(Matrix(mat,sparse = T))
  return(Matrix(mat))
}

#'
#'
avgDistancePoints<-function(points,varcol,dmat=NULL){
  show(paste0('Avg dist for ',nrow(points),' points'))
  if(nrow(points)<5000){sppoints = points}else{
    th = quantile(points[[varcol]],c(1 - 5000/nrow(points)))
    sppoints = points[points[[varcol]]>th,]
  }
  v = sppoints[[varcol]]
  if(is.null(dmat)){distances <- drop_units(st_distance(sppoints))}else{distances = dmat}
  #distances[distances==Inf]=0
  #distances[distances<0.01]=0
  #vv = (Diagonal(x=v/sum(v^2))%*%Matrix(distances,sparse=T))%*%Diagonal(x=v/sum(v^2))
  vv = (Diagonal(x=v)%*%Matrix(distances))%*%Diagonal(x=v)
  return(sum(vv)/(sum(v)^2*max(distances)))
}

entropyPoints<-function(points,varcol){
  v = points[[varcol]]
  v=v[v>0]
  v = v / sum(v)
  return(-sum(v*log(v))/log(length(v)))
}

slopePoints<-function(points,varcol){
  v = points[[varcol]]
  v=v[v>0]
  v=log(sort(v,decreasing=TRUE))
  rank = log(1:length(v))
  reg = lm(v~rank,data.frame(rank,v))
  return(reg$coefficients[2])
}


#########
## Raster functions


# weights for Moran
spatialWeights <- function (N,P,rescaleFactor=1){
  d_hor = (matrix(rep(cumsum(matrix(1,2*N+1,1)),2*P+1),nrow=2*N+1) - N - 1) ^ 2
  d_ver = (matrix(rep(cumsum(matrix(1,2*P+1,1)),2*N+1),ncol=2*P+1,byrow=TRUE) - P - 1) ^ 2
  w = 1 / (rescaleFactor*sqrt((d_hor^2 + d_ver^2)))
  w[w==Inf]=0
  return(w)
}



# total population
summaryPopulation<-function(m=NULL){
  if(is.null(m)){return(list(totalPop=NA,maxPop=NA,minPop=NA))}
  r_pop=raster(m)
  return(list(
    totalPop = cellStats(r_pop,'sum'),
    maxPop = cellStats(r_pop,'max'),
    minPop = cellStats(r_pop,'min')
  )
  )
}



#moran index 
moranIndex <- function(m=NULL){
  if(is.null(m)){return(list(moran=NA))}
  r_dens=raster(m/sum(m))
  return(list(moran=Moran(r_dens,spatialWeights(nrow(r_dens)-1,ncol(r_dens)-1))))
}

# same with use of focal
convolMoran <- function(r_pop){
  meanPop = cellStats(r_pop,sum)/ncell(r_pop)
  w = spatialWeights(nrow(r_pop)-1,ncol(r_pop)-1)
  return(ncell(r_pop) * cellStats(focal(r_pop-meanPop,w,sum,pad=TRUE,padValue=0)*(r_pop - meanPop),sum) / cellStats((r_pop - meanPop)*(r_pop - meanPop),sum) / cellStats(focal(raster(matrix(data=rep(1,ncell(r_pop)),nrow=nrow(r_pop))),w,sum,pad=TRUE,padValue=0),sum))
}


# average distance between indivduals
# normalized by max distance = world diagonal

distanceSubMatrix <- function(n_rows,n_cols,raster_cols,row_offset,col_offset){
  di = matrix(rep(cumsum(rep(1,n_rows)),n_cols),nrow=n_rows) + row_offset
  dj = matrix(rep(cumsum(rep(1,n_cols)),n_rows),nrow=n_rows,byrow=TRUE) + col_offset
  return(sqrt((abs(di-dj) %/%  raster_cols)^2 + (abs(di-dj) %% raster_cols)^2))
}

distanceMatrix <- function(N,P){
  d_hor = (matrix(rep(cumsum(matrix(1,2*N+1,1)),2*P+1),nrow=2*N+1) - N - 1) ^ 2
  d_ver = (matrix(rep(cumsum(matrix(1,2*P+1,1)),2*N+1),ncol=2*P+1,byrow=TRUE) - P - 1) ^ 2
  w = sqrt(d_hor + d_ver )
  return(w)
}


# average distance
# still very heavy computationally
# uses focal instead as in Moran Index computation.
#
averageDistance <- function(m=NULL){
  if(is.null(m)){return(list(averageDistance=NA))}
  r_pop=raster(m)
  return(list(averageDistance=cellStats(focal(r_pop,distanceMatrix(nrow(r_pop)-1,ncol(r_pop)-1),sum,pad=TRUE,padValue=0)*r_pop,sum) / ( cellStats(r_pop,sum)^2 * sqrt(nrow(r_pop)*ncol(r_pop)/pi))))
}


# distribution entropy --> rough equivalent of integrated local density ?
entropy <- function(m=NULL){
  if(is.null(m)){return(list(entropy=NA))}
  r_dens=raster(m/sum(m))
  m= values(r_dens)*cellStats(r_dens,function(x,...){na.omit(log(x))})
  m[is.na(m)]=0
  return(list(entropy = -1 / log(ncell(r_dens)) * sum(m) ))
}


# rank-size slope
# -> linear regression on sorted log series
rankSizeSlope <- function(m=NULL){
  if(is.null(m)){return(list(rankSizeAlpha =NA,rankSizeRSquared = NA))}
  r_pop = raster(m)
  size = cellStats(r_pop,function(x,...){na.omit(log(x))})
  size = size[size>0] # at least one person
  size=sort(size,decreasing=TRUE)
  #size = size[1:(length(size)*0.5)] # kill last quartile
  rank = log(1:length(size))
  if(length(size)>0){
    reg = lm(size~rank,data.frame(rank,size))
    return(list(
      rankSizeAlpha = reg$coefficients[2],
      rankSizeRSquared = summary(reg)$r.squared
    )
    )
  }else{return(list(rankSizeAlpha =NA,rankSizeRSquared = NA))}
}





# aggregate by resFactor (resolution = resolution * resFactor)
# a given square matrix of size areasize
simplifyBlock<-function(data,resFactor,areasize){
  m = matrix(data=data,nrow=areasize,byrow=TRUE)
  m[is.na(m)] <- 0
  res=matrix(0,areasize*resFactor,areasize*resFactor)
  for(x in 1:(areasize*resFactor)){
    for(y in 1:(areasize*resFactor)){
      res[x,y]=sum(m[((x-1)/resFactor+1):(x/resFactor),((y-1)/resFactor+1):(y/resFactor)])
    }
  }
  return(res)
}




# function to extract square subraster of a large raster
# used for visualisation of parts as real config
# - stored as temp asc file
extractSubRaster<- function(file,r,c,size,factor){
  raw <- raster(file)
  return(data.frame(simplifyBlock(getValuesBlock(raw,row=r,nrows=size,col=c,ncols=size),factor,size)/100))
  #r<-setExtent(r,extent(raster(paste0(Sys.getenv("CN_HOME"),'/Models/Synthetic/Density/temp_raster_pop.asc'))))
  #writeRaster(r,paste0(Sys.getenv("CN_HOME"),'/Models/Synthetic/Density/temp_raster_pop.asc'),format="ascii",overwrite=TRUE)
}

getCoordinates<-function(file,r,c){
  raw <- raster(file)
  return(spTransform(xyFromCell(raw,cellFromRowCol(raw,rownr = r,colnr = c),spatial = T),CRS("+proj=longlat +datum=WGS84")))
}



#system('rm -rf morphology')
#package.skeleton(name='morphology')

#install.packages(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Models/MorphologyPackage/morphology/'),repos=NULL,type='source')

