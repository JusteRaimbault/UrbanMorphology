
library(dplyr)
#library(sp)
#library(rgeos)
library(sf)
#library(tidyverse)
#library(ggplot2)
library(Matrix)

#source(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Models/MorphologyPackage/morphology.R'))
# for oml : must be in the same dir
source('morphology.R')

#library(devtools)
#devtools::install_github("tidyverse/ggplot2")
#require(ggplot2)

#'
#'
conditionalPercolation <- function(d,radius,popthq,nwcol,nwthq,gamma,decay,
                                   minclustsize=25,
                                   xcol="lonmin",ycol="latmin",popcol="totalPop",
                                   withMaps=T,
                                   resdir=paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Results/Percolation/Maps2/')
                                   ){
  
  if(withMaps==T){dir.create(resdir)}
  #popthq=0.95;nwthq=0.95;radius=50000;nwcol="euclPerf";popcol="totalPop";minclustsize=5
  #xcol="lonmin";ycol="latmin";popcol="totalPop";nwcol="mu";d=indics
  
  show(paste0("Parameter : radius = ",radius,' ; popthq = ',popthq,' ; nwcol = ',nwcol,' ; nwthq = ',nwthq,' ; gamma = ',gamma,' ; decay = ',decay))
  
  popth = quantile(d[[popcol]],c(popthq));nwth=quantile(d[[nwcol]],c(nwthq))
  
  show(paste0("Values : popth = ",popth));show(paste0("nwth = ",nwth));
  
  sppoints <- st_as_sf(d, coords = c(xcol, ycol), crs = 4326) %>% st_transform(3035)
  sppoints$index=1:nrow(sppoints)
  clusters = rep(NA,nrow(sppoints))
  currentclusternumber=1
  remainingpoints=T
  
  while(remainingpoints==T){
 
    cpoints = sppoints[is.na(clusters)&sppoints[[popcol]]>popth&sppoints[[nwcol]]>nwth,]
    show(nrow(cpoints))
    
    
    if(nrow(cpoints)>0){
      #' TODO this selection function should be passed as an argument -> can eg test random selection
      #' Rq : the process is independant of the choice of initial points
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
    if(nrow(currentcluster)>minclustsize){
      clusters[sppoints[currentcluster,op=st_intersects][["index"]]]=currentclusternumber
      currentclusternumber=currentclusternumber+1
    }else{
      clusters[sppoints[currentcluster,op=st_intersects][["index"]]]=-1
    }
    
    # this is wrong : could have a next cluster larger than minclustsize
    #remainingpoints=(nrow(cpoints)>0)&(nrow(currentcluster)>minclustsize)
    remainingpoints=(nrow(cpoints)>0)
    }else{
      remainingpoints=F
    }
    
  }
  
  clusters[clusters==-1]=NA
  sppoints$cluster = clusters
  
  #show(sppoints)
  #show(clusters)
  
  if(withMaps==T){
    #plot(sppoints%>%transmute(cluster=cluster))
    if(length(which(!is.na(clusters)))>0){
      g=ggplot(data.frame(x=indics$lonmin,y=indics$latmin,cluster=clusters),aes(x=x,y=y,fill=as.character(cluster)))
      ggsave(g+geom_raster()+scale_fill_discrete(guide=FALSE),file=paste0(resdir,popcol,popth,'_',nwcol,nwth,'_radius',radius,'.png'),width=15,height=10,units='cm')
    }
  }
    
  return(computeIndics(sppoints,gamma,decay))
  #return(sppoints)
}


#'
#' Aggregatd indices
aggregIndics <- function(res){
  show(paste0('Aggregating indicators for result with ',length(res$areas),' clusters'))
  return(list(
    area=sum(res$areas),
    pop=sum(res$areas),
    moran=mean(res$morans),
    avgdist=mean(res$avgdists),
    entropy=mean(res$entropies),
    slope=mean(res$slopes),
    efficiency=sum(efficiencies),
    emissions=sum(emissions)
  ))
}



#'
#'
computeIndics <- function(sppoints,gamma,decay,popcol="totalPop",emissioncol="CO2"){
  show(paste0('Computing indicators for ',nrow(sppoints),'points ; gamma = ',gamma,' ; decay = ',decay))
  clusteredpoints=sppoints[!is.na(sppoints$cluster),]
  areas=c();pops=c();morans=c();avgdists=c();entropies=c();slopes=c();
  efficiencies=c();emissions=c()
  clustnums = unique(clusteredpoints$cluster)
  for(i in 1:length(clustnums)){
    k=clustnums[i]
    show(paste0('Cluster ',i,'/',length(clustnums)))
    cluster = clusteredpoints[clusteredpoints$cluster==k,]
    show(paste0(length(cluster),' points'))
    envelope=st_convex_hull(st_union(cluster))
    #show(envelope)
    areas=append(areas,as.numeric(st_area(envelope)))
    allpoints=sppoints[envelope,op=function(x,y){st_is_within_distance(x,y,dist=10)}]
    
    #allpoints=sppoints[envelope,op=st_contains]
    currentpop=sum(allpoints[[popcol]]);show(paste0('pop = ',currentpop));pops=append(pops,currentpop)
    currentmoran = moranPoints(allpoints,popcol);show(paste0('moran = ',currentmoran));morans=append(morans,currentmoran)
    currentdist = avgDistancePoints(allpoints,popcol);show(paste0('avgdist = ',currentdist));avgdists=append(avgdists,currentdist)
    currententropy = entropyPoints(allpoints,popcol);show(paste0('entropy = ',currententropy));entropies=append(entropies,currententropy)
    currentslope = slopePoints(allpoints,popcol);show(paste0('slope = ',currentslope));slopes=append(slopes,currentslope)
    currentefficiency = sum(interactionPotential(allpoints,gravityGamma=gamma,gravityDistance=decay,weightcol=popcol));show(paste0('efficiency = ',currentefficiency))
    efficiencies=append(efficiencies,currentefficiency)
    currentemissions = sum(interactionPotential(allpoints,gravityGamma=gamma,gravityDistance=decay,weightcol=emissioncol));show(paste0('emissions = ',currentemissions))
    emissions=append(emissions,currentemissions)
  }
  return(list(areas=areas,pops=pops,morans=morans,avgdists=avgdists,
              entropies=entropies,slopes=slopes,efficiencies=efficiencies,emissions=emissions))
}



#' Compute a gravity potential between spatial points with population, of the form
#'   (P_i P_j / P_tot)^gamma * exp(-d_{ij} / d_0)
#'   TODO : influence of normalizing by pmax vs ptot ? -> investigate in interaction models
#'   - normalized such that sum = 1 ? NO to be able to have a size effects in comparing configurations
interactionPotential <- function(points,gravityGamma=1,gravityDistance=1,weightcol="totalPop"){
  show(paste0('Interaction potential for ',nrow(points),' points ; gravityGamma = ',gravityGamma,' ; gravityDistance = ',gravityDistance))
  # need a sampling if too much points -> take points with max value of weight
  if(nrow(points)<10000){sppoints = points}else{
    th = quantile(points[[weightcol]],c(1 - 10000/nrow(points)))
    sppoints = points[points[[weightcol]]>th,]
  }
  
  dmat = exp(-drop_units(st_distance(sppoints))/gravityDistance)
  diag(dmat)<-0;
  #dmat[dmat<0.01]=0
  #dmat = Matrix(dmat,sparse=T)
  #pops = sppoints[[weightcol]]
  #popmat = Diagonal(x=(pops/sum(pops))^gravityGamma)
  dmat = Matrix(dmat)
  weights=sppoints[[weightcol]]
  weights=weights/sum(weights) # ensure normalization
  # weights assumed already normalized
  weightmat=Diagonal(x=weights^gravityGamma)
  potentials = (weightmat%*%dmat)%*%weightmat
  return(as.matrix(potentials))
}









