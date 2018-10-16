
library(dplyr)
#library(sp)
#library(rgeos)
library(sf)
#library(tidyverse)
library(ggplot2)
library(Matrix)
library(igraph)


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
    # totally wrong assumption that cluster will be ordered by size..
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
#' Note : determinstic graph assumes deterministic percolation.
#'   precise assumptions of road nw percolation
graphPercolation <- function(d,radius,popthq,nwcol,nwthq,gamma,decay,
                             minclustsize=25,
                             xcol="lonmin",ycol="latmin",popcol="totalPop",
                             withMaps=T,
                             resdir=paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Results/Percolation/Maps2/')
){

  if(withMaps==T){dir.create(resdir)}
  show(paste0("Parameters : radius = ",radius,' ; popthq = ',popthq,' ; nwcol = ',nwcol,' ; nwthq = ',nwthq,' ; gamma = ',gamma,' ; decay = ',decay))
  popth = quantile(d[[popcol]],c(popthq));nwth=quantile(d[[nwcol]],c(nwthq))
  show(paste0("Values : popth = ",popth," ; nwth = ",nwth));
  sppoints <- st_as_sf(d, coords = c(xcol, ycol), crs = 4326) %>% st_transform(3035)
  
  cpoints = sppoints[sppoints[[popcol]]>popth&sppoints[[nwcol]]>nwth,]
  cpoints$cluster = rep(NA,nrow(cpoints))
  
  dmat = drop_units(st_distance(cpoints))
  adjmat = Matrix(dmat<=radius,sparse=T)
  g = graph_from_adjacency_matrix(adjmat,mode='undirected')
  comps = components(g)
  k=1
  for(cnum in which(comps$csize >= minclustsize)){
    cpoints$cluster[comps$membership==cnum]=k
    k=k+1
  }
  
  clusters=rep(NA,nrow(sppoints))
  clusters[sppoints[[popcol]]>popth&sppoints[[nwcol]]>nwth]=cpoints$cluster
  sppoints$cluster = clusters
  
  if(withMaps==T){
    #plot(sppoints%>%transmute(cluster=cluster))
    if(length(which(!is.na(clusters)))>0){
      g=ggplot(data.frame(x=d$lonmin,y=d$latmin,cluster=clusters),aes(x=x,y=y,fill=as.character(cluster)))
      ggsave(
        g+geom_raster()+scale_fill_discrete(guide=FALSE)+xlab("")+ylab("")+ggtitle(bquote(theta[P]*"="*.(popthq)*" ; "*theta[N]*"="*.(nwthq)*" ; "*.(nwcol)))
        ,file=paste0(resdir,popcol,popth,'_',nwcol,nwth,'_radius',radius,'.png'),width=20,height=18,units='cm')
    }
  }
  
  return(computeIndics(sppoints,gamma,decay))
  
}

#'
#' Aggregatd indices
aggregIndics <- function(res){
  #show(paste0('Aggregating indicators for result with ',length(res$areas),' clusters'))
  return(c(
    area=sum(res$areas),
    pop=sum(res$pops),
    avgpop=mean(res$pops),
    pophierarchy=hierarchy(res$pops),
    moran=mean(res$morans),
    avgdist=mean(res$avgdists),
    entropy=mean(res$entropies),
    slope=mean(res$slopes),
    efficiency=sum(res$efficiencies),
    avgefficiency=mean(res$efficiencies),
    emissions=sum(res$emissions),
    avgemissions=mean(res$emissions),
    totalemissions=sum(res$totalemissions),
    avgtotemissions=mean(res$totalemissions)
  ))
}

hierarchy <- function(x){
  x=x[x>0]
  if(length(x)>1){
    #show(log(1:length(x)));show(log(sort(x,decreasing = T)))
    d=data.frame(rank=log(1:length(x)),size=log(sort(x,decreasing = T)))
  reg = lm(size~rank,data=d)
  return(reg$coefficients[2])
  }else{return(0)}
}


#'
#'
computeIndics <- function(sppoints,gamma,decay,dmat=NULL,popcol="totalPop",emissioncol="CO2"){
  show(paste0('Computing indicators for ',nrow(sppoints),'points ; gamma = ',gamma,' ; decay = ',decay))
  #clusteredpoints=sppoints[!is.na(sppoints$cluster),]
  areas=c();pops=c();morans=c();avgdists=c();entropies=c();slopes=c();
  efficiencies=c();emissions=c();totalemissions=c()
  clustnums = unique(sppoints$cluster);clustnums=clustnums[!is.na(clustnums)]
  if(length(clustnums)>0){
  for(i in 1:length(clustnums)){
    k=clustnums[i]
    show(paste0('Cluster ',i,'/',length(clustnums)))
    #cluster = clusteredpoints[clusteredpoints$cluster==k,]
    inds = which(sppoints$cluster==k)
    cluster = sppoints[inds,]
    show(paste0(nrow(cluster),' points'))
    
    # do not include convex area ? can be useful for high nw ths
    envelope=st_convex_hull(st_union(cluster))
    #show(envelope)
    areas=append(areas,as.numeric(st_area(envelope)))
    #allpoints=sppoints[envelope,op=function(x,y){st_is_within_distance(x,y,dist=10)}]
    allpoints=sppoints[envelope,op=st_within]
    
    if(nrow(allpoints)>0){
      #allpoints=sppoints[envelope,op=st_contains]
      currentpop=sum(allpoints[[popcol]]);show(paste0('pop = ',currentpop));pops=append(pops,currentpop)
      currentmoran = moranPoints(allpoints,popcol);show(paste0('moran = ',currentmoran));morans=append(morans,currentmoran)
      currentdist = avgDistancePoints(allpoints,popcol);show(paste0('avgdist = ',currentdist));avgdists=append(avgdists,currentdist)
      currententropy = entropyPoints(allpoints,popcol);show(paste0('entropy = ',currententropy));entropies=append(entropies,currententropy)
      currentslope = slopePoints(allpoints,popcol);show(paste0('slope = ',currentslope));slopes=append(slopes,currentslope)
      currentefficiency = sum(interactionPotential(allpoints,gravityGamma=gamma,gravityDistance=decay,weightcol=popcol));show(paste0('efficiency = ',currentefficiency))
      efficiencies=append(efficiencies,currentefficiency)
      currentemissions = sum(interactionPotential(allpoints,gravityGamma=gamma,gravityDistance=decay,weightcol=emissioncol));show(paste0('emissions = ',currentemissions))
      emissions=append(emissions,currentemissions);
      totalemissions=append(totalemissions,sum(allpoints[[emissioncol]]))
    }
  }
  }else{
    areas=c(NA);pops=c(NA);morans=c(NA);avgdists=c(NA);entropies=c(NA);slopes=c(NA);
    efficiencies=c(NA);emissions=c(NA);totalemissions=c(NA)
  }
  return(list(areas=areas,pops=pops,morans=morans,avgdists=avgdists,
              entropies=entropies,slopes=slopes,efficiencies=efficiencies,emissions=emissions,totalemissions=totalemissions))
}



#' Compute a gravity potential between spatial points with population, of the form
#'   (P_i P_j / P_tot)^gamma * exp(-d_{ij} / d_0)
#'   TODO : influence of normalizing by pmax vs ptot ? -> investigate in interaction models
#'   - normalized such that sum = 1 ? NO to be able to have a size effects in comparing configurations
interactionPotential <- function(points,gravityGamma=1,gravityDistance=1,weightcol="totalPop",dmat=NULL){
  show(paste0('Interaction potential for ',nrow(points),' points ; gravityGamma = ',gravityGamma,' ; gravityDistance = ',gravityDistance))
  # need a sampling if too much points -> take points with max value of weight
  if(nrow(points)<5000){sppoints = points}else{
    th = quantile(points[[weightcol]],c(1 - 5000/nrow(points)))
    sppoints = points[points[[weightcol]]>th,]
  }
  
  if(is.null(dmat)){dmat = drop_units(st_distance(sppoints))}
  
  wmat = exp(-dmat/gravityDistance)
  diag(wmat)<-0;
  #dmat[dmat<0.01]=0
  #dmat = Matrix(dmat,sparse=T)
  #pops = sppoints[[weightcol]]
  #popmat = Diagonal(x=(pops/sum(pops))^gravityGamma)
  wmat = Matrix(wmat)
  weights=sppoints[[weightcol]]
  weights=weights/sum(weights) # ensure normalization
  # weights assumed already normalized
  weightmat=Diagonal(x=weights^gravityGamma)
  potentials = (weightmat%*%wmat)%*%weightmat
  return(as.matrix(potentials))
}









