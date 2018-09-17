
library(dplyr)
#library(sp)
#library(rgeos)
library(sf)
library(tidyverse)
library(ggplot2)
library(Matrix)

#library(devtools)
#devtools::install_github("tidyverse/ggplot2")
#require(ggplot2)

#'
#'
conditionalPercolation <- function(d,popthq,nwthq,radius,minclustsize=5,
                                   xcol="lonmin",ycol="latmin",popcol="totalPop",nwcol="euclPerf",
                                   resdir=paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Results/Percolation/Maps/')
                                   ){
  dir.create(resdir)
  #popth=5000000;nwth=3000;radius=50000
  #xcol="lonmin";ycol="latmin";popcol="totalPop";nwcol="mu";d=indics
  #x=data[,xcol];y=data[,ycol];population=data[,popcol];network=data[,nwcol]
  #sppoints = SpatialPoints(data.frame(lon=x,lat==y),proj4string = countries@proj4string)
  
  show(paste0("radius = ",radius))
  
  popth = quantile(d[[popcol]],c(popthq));nwth=quantile(d[[nwcol]],c(nwthq))
  
  show(paste0("popth = ",popth));show(paste0("nwth = ",nwth));
  
  sppoints <- st_as_sf(d, coords = c(xcol, ycol), crs = 4326) %>% st_transform(3035)
  sppoints$index=1:nrow(sppoints)
  clusters = rep(NA,nrow(sppoints))
  currentclusternumber=1
  remainingpoints=T
  
  while(remainingpoints==T){
 
    cpoints = sppoints[is.na(clusters)&sppoints[[popcol]]>popth&sppoints[[nwcol]]>nwth,]
    show(nrow(cpoints))
    
    
    if(nrow(cpoints)>0){
      currentcluster = cpoints[cpoints[[popcol]]==max(cpoints[[popcol]]),]
    newpoints=T
    while(newpoints==T){
      buffer = st_buffer(currentcluster, dist = radius)
      prevpoints = nrow(currentcluster)
      currentcluster <- cpoints[buffer,op=st_within]
      newpoints=(nrow(currentcluster) - prevpoints > 0)
      #show(nrow(currentcluster))
    }
    #plot(currentcluster,max.plot=1)
    # get indices in all points
    clusters[sppoints[currentcluster,op=st_intersects][["index"]]]=currentclusternumber
    currentclusternumber=currentclusternumber+1
    
    remainingpoints=(nrow(cpoints)>0)&(nrow(currentcluster)>minclustsize)
    }else{
      remainingpoints=F
    }
    
  }
  
  sppoints$cluster = clusters
  
  #show(sppoints)
  #show(clusters)
  
  #plot(sppoints%>%transmute(cluster=cluster))
  if(length(which(!is.na(clusters)))>0){
    g=ggplot(data.frame(x=indics$lonmin,y=indics$latmin,cluster=clusters),aes(x=x,y=y,fill=as.character(cluster)))
    ggsave(g+geom_raster()+scale_fill_discrete(name="cluster"),file=paste0(resdir,popcol,popth,'_',nwcol,nwth,'_radius',radius,'.png'),width=15,height=10,units='cm')
  }
  
  #return(computeClustersIndics(sppoints))
  return(sppoints)
}


#'
#'
computeIndics <- function(sppoints,popcol="totalPop"){
  
}




#'
#'
computeClustersIndics <- function(sppoints,popcol="totalPop"){
  clusteredpoints=sppoints[!is.na(sppoints$cluster),]
  areas=c();pops=c()
  for(k in unique(clusteredpoints$cluster)){
    show(k)
    cluster = clusteredpoints[clusteredpoints$cluster==k,]
    envelope=st_convex_hull(st_union(cluster))
    areas=append(areas,as.numeric(st_area(envelope)))
    allpoints=sppoints[envelope,op=function(x,y){st_is_within_distance(x,y,dist=10)}]
    #allpoints=sppoints[envelope,op=st_contains]
    pops=append(pops,sum(allpoints[[popcol]]))
  }
  return(list(areas=areas,pops=pops))
}


#'
#' Compute morphological indices
morphoIndices <-function(sppoints,popcol="totalPop"){
  
}
  


#' Compute a gravity potential between spatial points with population, of the form
#'   (P_i P_j / P_tot)^gamma * exp(-d_{ij} / d_0)
#'   TODO : influence of normalizing by pmax vs ptot ? -> investigate in interaction models
#'   
#'   - normalized such that sum = 1
interactionPotential <- function(sppoints,gravityGamma=1,gravityDistance=1,popcol="totalPop"){
  dmat = exp(-st_distance(sppoints)/gravityDistance)
  diag(dmat)<-0
  dmat = Matrix(dmat)
  pops = sppoints[[popcol]]
  popmat = Diagonal(x=(pops/sum(pops))^gravityGamma)
  potentials = popmat%*%dmat%*%popmat
  return(as.matrix(potentials))
}









