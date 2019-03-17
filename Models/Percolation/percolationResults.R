
setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Models/Percolation'))

# indics
load('../../Data/consolidated/indics.RData')

# load results from percolation

#resstamp = '20181016_024956'
#purpose='directSampling'
#resstamp = '20190316_182550'
resstamp='20190316_235006'
purpose='directSampling_radius'

#resstamp = 'test_20190316_154614'
load(paste0('res/',purpose,'_',resstamp,'.RData'))
load(paste0('res/',purpose,'_params_',resstamp,'.RData'))
resdir = paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology/Results/Percolation/',purpose,'_',resstamp,'/');dir.create(resdir)


source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))


totalPop = sum(indics$totalPop)

#####
## full areas stats

for(k in 1:length(res)){res[[k]]$paramrow=rep(k,length(res[[k]]$pops))}
full = data.frame(matrix(unlist(lapply(res,function(l){l$areas=NULL;l$clustdiameters=NULL;l$clustsizesnodes=NULL;l$clustsizesedges=NULL;unlist(data.frame(t(data.frame(l))))})),byrow=T,ncol=length(res[[1]])-4))
colnames(full)<-names(res[[1]])[c(-1,-10,-11,-12)]
params$paramrow=1:nrow(params)
full=left_join(full,params)

# compute relative indicators
full$relemissions=rep(nrow(full),0)
full$relefficiency=rep(nrow(full),0)
for(k in unique(full$paramrow)){
  full$relemissions[full$paramrow==k]=(max(full$emissions[full$paramrow==k])-full$emissions[full$paramrow==k])/(max(full$emissions[full$paramrow==k])-min(full$emissions[full$paramrow==k]))
  full$relefficiency[full$paramrow==k]= 1 - (max(full$efficiencies[full$paramrow==k])-full$efficiencies[full$paramrow==k])/(max(full$efficiencies[full$paramrow==k])-min(full$efficiencies[full$paramrow==k])) 
}


#####
# percolation transition

#sapply(res[[1]],length)
#sapply(res,length)
#sapply(res,function(l){length(l[[1]])})
#max(sapply(res,function(l){abs(length(l[[1]])-length(l[[10]]))+abs(length(l[[10]])-length(l[[11]]))+abs(length(l[[11]])-length(l[[12]]))+abs(length(l[[12]])-length(l[[13]]))}))
# which(params$nwcol=='mu'&params$nwthq==0.95&params$radius==75000)


# -> need to filter NAs
# BEFORE the lapply (not same length) - ULTRA DIRTY data structure -> really need to stop using R...
resfilt = list();paraminds=c();kk=1;for(k in 1:length(res)){if(length(which(is.na(res[[k]][[1]])))==0){resfilt[[kk]]=res[[k]];paraminds=append(paraminds,k);kk=kk+1}}
res=resfilt
params=params[paraminds,]

toremove=c("pops","morans","avgdists","entropies","slopes","efficiencies","emissions","totalemissions")
for(k in 1:length(res)){res[[k]]$paramrow=rep(k,length(res[[k]]$areas))}
fulltrans = data.frame(matrix(unlist(lapply(res,function(l){for(tr in toremove){l[[tr]]=NULL};unlist(data.frame(t(data.frame(l))))})),byrow=T,ncol=length(res[[1]])-length(toremove)))
colnames(fulltrans)<-names(res[[1]])[c(-9:-2)]
params$paramrow=1:nrow(params)
fulltrans=as.tbl(left_join(fulltrans,params))

dir.create(paste0(resdir,'transition'))

# not normalized
fulltrans$nwparam = paste0(fulltrans$nwcol,fulltrans$nwthq)
# maximum within a parameter value
maxs = as.tbl(fulltrans)%>% group_by(nwcol,popthq,nwthq,radius,nwparam)%>%summarize(maxsizenodes = max(clustsizesnodes),maxsizeedges=max(clustsizesedges))

g=ggplot(maxs[maxs$nwthq>0|(maxs$nwthq==0&maxs$nwcol=='ecount'),],aes(x=radius,y=maxsizenodes,group=nwparam,color=nwparam))
g+geom_point()+geom_line()+xlab(expression(r[0]*' (m)'))+ylab("Size of largest cluster (nodes)")+
  scale_color_discrete(name='Parameters',labels=c(expression(theta[N]*'= 0'),expression(theta[N]*'= 0.8 ; ecount'),expression(theta[N]*'= 0.8 ; euclPerf'),expression(theta[N]*'= 0.95 ; euclPerf'),expression(theta[N]*'= 0.8 ; mu'),expression(theta[N]*'= 0.95 ; mu'),expression(theta[N]*'= 0.8 ; vcount'),expression(theta[N]*'= 0.95 ; vcount')))+
  stdtheme
ggsave(file=paste0(resdir,'transition/abssize_nodes.png'),width=30,height=20,units='cm')


g=ggplot(maxs[maxs$nwthq>0|(maxs$nwthq==0&maxs$nwcol=='ecount'),],aes(x=radius,y=maxsizeedges,group=nwparam,color=nwparam))
g+geom_point()+geom_line()+xlab(expression(r[0]*' (m)'))+ylab("Size of largest cluster (edges)")+
  scale_color_discrete(name='Parameters',labels=c(expression(theta[N]*'= 0'),expression(theta[N]*'= 0.8 ; ecount'),expression(theta[N]*'= 0.8 ; euclPerf'),expression(theta[N]*'= 0.95 ; euclPerf'),expression(theta[N]*'= 0.8 ; mu'),expression(theta[N]*'= 0.95 ; mu'),expression(theta[N]*'= 0.8 ; vcount'),expression(theta[N]*'= 0.95 ; vcount')))+
  stdtheme
ggsave(file=paste0(resdir,'transition/abssize_edges.png'),width=30,height=20,units='cm')


# overall maximums
maxall = as.tbl(fulltrans)%>% group_by(nwcol,popthq,nwthq,nwparam)%>%summarize(maxsizenodes = max(clustsizesnodes),maxsizeedges=max(clustsizesedges))

# local maximums
maxs = as.tbl(fulltrans)%>% group_by(nwcol,popthq,nwthq,radius,nwparam)%>%summarize(relsizenodes = max(clustsizesnodes)/maxall$maxsizenodes[maxall$nwparam==nwparam[1]],relsizeedges=max(clustsizesedges)/maxall$maxsizeedges[maxall$nwparam==nwparam[1]])

g=ggplot(maxs[maxs$nwthq>0|(maxs$nwthq==0&maxs$nwcol=='ecount'),],aes(x=radius,y=relsizenodes,group=nwparam,color=nwparam))
g+geom_point()+geom_line()+xlab(expression(r[0]*' (m)'))+ylab("Relative size of largest cluster (nodes)")+
  scale_color_discrete(name='Parameters',labels=c(expression(theta[N]*'= 0'),expression(theta[N]*'= 0.8 ; ecount'),expression(theta[N]*'= 0.8 ; euclPerf'),expression(theta[N]*'= 0.95 ; euclPerf'),expression(theta[N]*'= 0.8 ; mu'),expression(theta[N]*'= 0.95 ; mu'),expression(theta[N]*'= 0.8 ; vcount'),expression(theta[N]*'= 0.95 ; vcount')))+
  stdtheme
ggsave(file=paste0(resdir,'transition/relsize_nodes.png'),width=30,height=20,units='cm')


g=ggplot(maxs,aes(x=radius,y=relsizeedges,group=nwparam,color=nwparam))
g+geom_point()+geom_line()+xlab(expression(r[0]*' (m)'))+ylab("Relative size of largest cluster (edges)")+
  scale_color_discrete(name='Parameters',labels=c(expression(theta[N]*'= 0'),expression(theta[N]*'= 0.8 ; ecount'),expression(theta[N]*'= 0.8 ; euclPerf'),expression(theta[N]*'= 0.95 ; euclPerf'),expression(theta[N]*'= 0.8 ; mu'),expression(theta[N]*'= 0.95 ; mu'),expression(theta[N]*'= 0.8 ; vcount'),expression(theta[N]*'= 0.95 ; vcount')))+
  stdtheme
ggsave(file=paste0(resdir,'transition/relsize_edges.png'),width=30,height=20,units='cm')


## fractal dimension
# estimated by log(# Nodes) = lm(log(diameter))

#counts = as.tbl(fulltrans)%>% group_by(nwcol,popthq,nwthq,radius,nwparam)%>%summarize(count=n())
# length(which(counts$count<10))
# d=fulltrans[fulltrans$nwcol=='mu'&fulltrans$nwthq==0.95&fulltrans$radius==75000,]
# reg = lm(data = data.frame(lnodes=log(d$clustsizesnodes),ldiam=log(d$clustdiameters)),lnodes~ldiam)

minobs=20

fractdim = as.tbl(fulltrans)%>% group_by(nwcol,popthq,nwthq,radius,nwparam)%>%summarize(
  fractdim=ifelse(n()<minobs,NA,summary(lm(data = data.frame(lnodes=log(clustsizesnodes),ldiam=log(clustdiameters)),lnodes~ldiam))$coefficients[2,1]),
  fractdimsd=ifelse(n()<minobs,NA,summary(lm(data = data.frame(lnodes=log(clustsizesnodes),ldiam=log(clustdiameters)),lnodes~ldiam))$coefficients[2,2]),
  fractdimadjrsquared=ifelse(n()<minobs,NA,summary(lm(data = data.frame(lnodes=log(clustsizesnodes),ldiam=log(clustdiameters)),lnodes~ldiam))$adj.r.squared),
  fractrelsd=fractdimsd/fractdim
)

g=ggplot(fractdim[(fractdim$nwthq>0|(fractdim$nwthq==0&fractdim$nwcol=='ecount'))&fractdim$radius<=50000,],aes(x=radius,y=fractdim,ymin=fractdim-fractdimsd,ymax=fractdim+fractdimsd,group=nwparam,color=nwparam))
g+geom_point()+#geom_errorbar()+
  geom_line()+xlab(expression(r[0]*' (m)'))+ylab("Fractal dimension")+
  scale_color_discrete(name='Parameters',labels=c(expression(theta[N]*'= 0'),expression(theta[N]*'= 0.8 ; ecount'),expression(theta[N]*'= 0.8 ; euclPerf'),expression(theta[N]*'= 0.95 ; euclPerf'),expression(theta[N]*'= 0.8 ; mu'),expression(theta[N]*'= 0.95 ; mu'),expression(theta[N]*'= 0.8 ; vcount'),expression(theta[N]*'= 0.95 ; vcount')))+
  stdtheme
ggsave(file=paste0(resdir,'transition/fractaldimension.png'),width=30,height=20,units='cm')
# -> shitty as work with geo diams ?

g+geom_point()+geom_errorbar()+
  geom_line()+xlab(expression(r[0]*' (m)'))+ylab("Fractal dimension")+
  scale_color_discrete(name='Parameters',labels=c(expression(theta[N]*'= 0'),expression(theta[N]*'= 0.8 ; ecount'),expression(theta[N]*'= 0.8 ; euclPerf'),expression(theta[N]*'= 0.95 ; euclPerf'),expression(theta[N]*'= 0.8 ; mu'),expression(theta[N]*'= 0.95 ; mu'),expression(theta[N]*'= 0.8 ; vcount'),expression(theta[N]*'= 0.95 ; vcount')))+
  stdtheme
ggsave(file=paste0(resdir,'transition/fractaldimension_errorbars.png'),width=30,height=20,units='cm')

summary(fractdim)
# signif of variations (min/max)
fractdimsignif = fractdim%>%filter(!is.na(fractdim))%>%group_by(nwcol,popthq,nwthq,nwparam)%>%summarize(
  signif = (max(fractdim)-fractdimsd[which(fractdim==max(fractdim))])-(min(fractdim)+fractdimsd[fractdim==min(fractdim)])
)


########
## Optimization on sustainibility indicators


# within config pareto fronts ?
popthq=c(0.85,0.95);nwthq=c(0.0,0.95);radius=c(8000,50000)
g=ggplot(full[full$popthq%in%popthq&full$nwthq%in%nwthq&full$radius%in%radius,],aes(x=relemissions,y=relefficiency,color=paste0(popthq," ; ",nwthq," ; ",radius),size= pops/totalPop))
g+geom_point(alpha=0.5)+facet_wrap(~paste0(decay," ; ",gamma),scales="free")+
  scale_color_discrete(name=expression(theta[P]*" ; "*theta[N]*" ; "*r[0]))+scale_size_continuous(name="Relative population")+
  xlab("Relative potential emissions")+ylab("1 - Relative potential efficiency")+stdtheme
ggsave(file=paste0(resdir,'full_paretos.png'),width=45,height=30,units='cm')


# link to morphology
morph=full[!is.na(full$morans),c("morans","avgdists","entropies","slopes")]
for(j in 1:ncol(morph)){morph[,j]=(max(morph[,j])-morph[,j])/(max(morph[,j])-min(morph[,j]))}
pca = prcomp(morph[!duplicated(morph),])
# Cumulative Proportion  0.7903 0.9570 0.98717 1.00000
#                   PC1        PC2         PC3        PC4
#morans    -0.05350415  0.0291656  0.22191261  0.9731606
#avgdists   0.69048222 -0.7190648 -0.03881322  0.0683636
#entropies  0.54025216  0.4684572  0.68479020 -0.1404914
#slopes     0.47801590  0.5124870 -0.69304452  0.1689589
# ~ level of monocentricity ?
full$PC1 = rep(NA,nrow(full));full$PC1[!is.na(full$morans)]=(as.matrix(morph)%*%pca$rotation)[,1]

popthq=c(0.85,0.95);nwthq=c(0.0,0.95);radius=c(8000,50000)
g=ggplot(full[full$popthq%in%popthq&full$nwthq%in%nwthq&full$radius%in%radius,],aes(x=PC1,y=relemissions,color=paste0(popthq," ; ",nwthq," ; ",radius)))
g+geom_point(alpha=0.5,pch=".")+geom_smooth(span=0.8)+facet_wrap(~paste0(decay," ; ",gamma),scales="free")+
  scale_color_discrete(name=expression(theta[P]*" ; "*theta[N]*" ; "*r[0]))+
  xlab("PC1")+ylab("Relative potential emissions")+stdtheme
ggsave(file=paste0(resdir,'full_pc1-emissions.png'),width=45,height=30,units='cm')

g=ggplot(full[full$popthq%in%popthq&full$nwthq%in%nwthq&full$radius%in%radius,],aes(x=PC1,y=relefficiency,color=paste0(popthq," ; ",nwthq," ; ",radius)))
g+geom_point(alpha=0.5,pch=".")+geom_smooth(span=0.8)+facet_wrap(~paste0(decay," ; ",gamma),scales="free")+
  scale_color_discrete(name=expression(theta[P]*" ; "*theta[N]*" ; "*r[0]))+
  xlab("PC1")+ylab("1 - Relative potential efficiency")+stdtheme
ggsave(file=paste0(resdir,'full_pc1-efficiency.png'),width=45,height=30,units='cm')


## paretos pop-emissions without gravity potentials
popthq=c(0.85,0.95);nwthq=c(0.0,0.95);radius=c(8000,50000)
g=ggplot(full[full$popthq%in%popthq&full$nwthq%in%nwthq&full$radius%in%radius&full$gamma==1&full$decay==100,],aes(x=log10(totalemissions),y=1-log10(pops),color=paste0(popthq," ; ",nwthq," ; ",radius)))
g+geom_point(alpha=0.5)+scale_color_discrete(name=expression(theta[P]*" ; "*theta[N]*" ; "*r[0]))+
  xlab("log(Total emissions)")+ylab("1 - log(Total population)")+stdtheme
ggsave(file=paste0(resdir,'full_effective_pareto.png'),width=20,height=15,units='cm')


# link with morpho
g=ggplot(full[full$popthq%in%popthq&full$nwthq%in%nwthq&full$radius%in%radius&full$gamma==1&full$decay==100,],aes(x=PC1,y=log10(pops),color=paste0(popthq," ; ",nwthq," ; ",radius)))
g+geom_point(pch='.')+geom_smooth()+scale_color_discrete(name=expression(theta[P]*" ; "*theta[N]*" ; "*r[0]))+
  xlab("PC1")+ylab("log(Total population)")+stdtheme
ggsave(file=paste0(resdir,'full_effective_morpho-pop.png'),width=20,height=15,units='cm')

g=ggplot(full[full$popthq%in%popthq&full$nwthq%in%nwthq&full$radius%in%radius&full$gamma==1&full$decay==100,],aes(x=PC1,y=log10(totalemissions),color=paste0(popthq," ; ",nwthq," ; ",radius)))
g+geom_point(pch='.')+geom_smooth()+scale_color_discrete(name=expression(theta[P]*" ; "*theta[N]*" ; "*r[0]))+
  xlab("PC1")+ylab("log(Total emissions)")+stdtheme
ggsave(file=paste0(resdir,'full_effective_morpho-emissions.png'),width=20,height=15,units='cm')


######
## aggregated stats

aggregated = data.frame(matrix(unlist(lapply(res,aggregIndics)),byrow = T,nrow=length(res)))
colnames(aggregated)<-names(aggregIndics(res[[1]]))
sres = cbind(params,aggregated)
sres=sres[!is.na(sres$pop),]

#g=ggplot(sres,aes(x=efficiency,y=emissions,col=pop/totalPop))
#g+geom_point()+facet_wrap(~paste0(decay," ; ",gamma),scales="free")
# -> not useful : look at normalized plot

# normalize for each gravity param
for(gamma in gammas){for (decay in decays){
  sres[sres$gamma==gamma&sres$decay==decay,"emissions"] = (max(sres[sres$gamma==gamma&sres$decay==decay,"emissions"]) - sres[sres$gamma==gamma&sres$decay==decay,"emissions"])/ (max(sres[sres$gamma==gamma&sres$decay==decay,"emissions"]) - min(sres[sres$gamma==gamma&sres$decay==decay,"emissions"]))
  sres[sres$gamma==gamma&sres$decay==decay,"efficiency"] =  1 - (max(sres[sres$gamma==gamma&sres$decay==decay,"efficiency"]) - sres[sres$gamma==gamma&sres$decay==decay,"efficiency"])/ (max(sres[sres$gamma==gamma&sres$decay==decay,"efficiency"]) - min(sres[sres$gamma==gamma&sres$decay==decay,"efficiency"]))
}}


g=ggplot(sres,aes(x=emissions,y=efficiency,col=pop/totalPop))
g+geom_point()+facet_grid(decay~gamma)+
  xlab("Normalized potential emissions")+ylab("1 - Normalized potential efficiency")+
  scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_paretos.png'),width=30,height=25,units='cm')

g=ggplot(sres,aes(x=emissions*totalPop/pop,y=efficiency*totalPop/pop,col=pop/totalPop))
g+geom_point()+facet_grid(decay~gamma)+xlab("Relative potential emissions")+ylab("Relative potential unefficiency")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_paretos_relative.png'),width=30,height=25,units='cm')



# test : which params produce pareto fronts ?
g=ggplot(sres,aes(x=emissions,y=efficiency,col=nwcol))
g+geom_point(alpha=0.5)+facet_grid(decay~gamma)+xlab("Normalized potential emissions")+ylab("1 - Normalized potential efficiency")+scale_color_discrete(name="Percolation\nIndicator")+stdtheme
ggsave(file=paste0(resdir,'aggreg_paretos_nwcol.png'),width=30,height=25,units='cm')

g=ggplot(sres,aes(x=emissions*totalPop/pop,y=efficiency*totalPop/pop,col=nwcol))
g+geom_point(alpha=0.5)+facet_grid(decay~gamma)+xlab("Relative potential emissions")+ylab("Relative potential unefficiency")+scale_color_discrete(name="Percolation\nIndicator")+stdtheme
ggsave(file=paste0(resdir,'aggreg_paretos_nwcol_relative.png'),width=30,height=25,units='cm')


g=ggplot(sres,aes(x=emissions,y=efficiency,col=radius,size=popthq))
g+geom_point(alpha=0.5)+facet_grid(decay~gamma)+xlab("Normalized potential emissions")+ylab("1 - Normalized potential efficiency")+scale_color_continuous(name=expression(r[0]))+scale_size_continuous(name=expression(theta[P]))+stdtheme
ggsave(file=paste0(resdir,'aggreg_paretos_radiuspopthq.png'),width=30,height=25,units='cm')

g=ggplot(sres,aes(x=emissions*totalPop/pop,y=efficiency*totalPop/pop,col=radius,size=popthq))
g+geom_point(alpha=0.5)+facet_grid(decay~gamma)+xlab("Relative potential emissions")+ylab("Relative potential unefficiency")+scale_color_continuous(name=expression(r[0]))+scale_size_continuous(name=expression(theta[P]))+stdtheme
ggsave(file=paste0(resdir,'aggreg_paretos_radiuspopthq_relative.png'),width=30,height=25,units='cm')


g=ggplot(sres,aes(x=emissions,y=efficiency,col=radius,size=nwthq))
g+geom_point(alpha=0.5)+facet_grid(decay~gamma)+xlab("Normalized potential emissions")+ylab("1 - Normalized potential efficiency")+scale_color_continuous(name=expression(r[0]))+scale_size_continuous(name=expression(theta[N]))+stdtheme
ggsave(file=paste0(resdir,'aggreg_paretos_radiusnwthq.png'),width=30,height=25,units='cm')

g=ggplot(sres,aes(x=emissions*totalPop/pop,y=efficiency*totalPop/pop,col=radius,size=nwthq))
g+geom_point(alpha=0.5)+facet_grid(decay~gamma)+xlab("Relative potential emissions")+ylab("Relative potential unefficiency")+scale_color_continuous(name=expression(r[0]))+scale_size_continuous(name=expression(theta[N]))+stdtheme
ggsave(file=paste0(resdir,'aggreg_paretos_radiusnwthq_relative.png'),width=30,height=25,units='cm')



g=ggplot(sres,aes(x=emissions,y=efficiency,col=as.character(nwthq),size=popthq))
g+geom_point(alpha=0.5)+facet_grid(decay~gamma)+xlab("Normalized potential emissions")+ylab("1 - Normalized potential efficiency")+scale_color_discrete(name=expression(theta[N]))+scale_size_continuous(name=expression(theta[P]))+stdtheme
ggsave(file=paste0(resdir,'aggreg_paretos_popthqnwthq.png'),width=30,height=25,units='cm')

g=ggplot(sres,aes(x=emissions*totalPop/pop,y=efficiency*totalPop/pop,col=as.character(nwthq),size=popthq))
g+geom_point(alpha=0.5)+facet_grid(decay~gamma)+xlab("Relative potential emissions")+ylab("Relative potential unefficiency")+scale_color_discrete(name=expression(theta[N]))+scale_size_continuous(name=expression(theta[P]))+stdtheme
ggsave(file=paste0(resdir,'aggreg_paretos_popthqnwthq_relative.png'),width=30,height=25,units='cm')



alpha = lm(y~x,data=data.frame(x=log(sres$pop/totalPop),y=log(sres$totalemissions)))$coefficients[2]
g=ggplot(sres,aes(x=pop/totalPop,y=totalemissions))
g+geom_point()+scale_x_log10()+scale_y_log10()+geom_smooth(method="lm")+
  ggtitle(bquote(alpha*"="*.(alpha)))+xlab("Relative population")+ylab("Total emissions")+stdtheme
ggsave(file=paste0(resdir,'aggreg_effective.png'),width=20,height=18,units='cm')



###
# cluster balance ? (territorial equity ?)
# : alpha pop already good measure.

# -> work also with cluster level measures, or other statistics.

# relative efficiency ?

# link morpho - perfs
#g=ggplot(sres,aes(x=moran,y=totalemissions,col=pop/totpop))
#g+geom_point()+geom_smooth(method='loess',n=50)
#g=ggplot(sres,aes(x=moran,y=emissions,col=pop/totpop))
#g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',n=50)

morph=sres[,c("moran","avgdist","entropy","slope")]
for(j in 1:ncol(morph)){morph[,j]=(max(morph[,j])-morph[,j])/(max(morph[,j])-min(morph[,j]))}
pca = prcomp(morph[!duplicated(morph),])
# summary(pca)
# Cumulative Proportion  0.7296 0.9650 0.99761 1.00000
#                PC1       PC2         PC3         PC4
#moran   -0.3088585 0.9493848 -0.04444327  0.03605266
#avgdist  0.5417362 0.1415668 -0.82239570 -0.10072759
#entropy  0.5108424 0.2140447  0.45847647 -0.69499942
#slope    0.5917502 0.1811415  0.33390034  0.71100630

sres$PC1 = (as.matrix(morph)%*%pca$rotation)[,1]
sres$PC2 = (as.matrix(morph)%*%pca$rotation)[,2]

g=ggplot(sres,aes(x=PC1,y=emissions,col=pop/totpop))
g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',n=50)+xlab("PC1")+ylab("Normalized potential emissions")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_pc1-emissions.png'),width=30,height=25,units='cm')

# targeted plot
g=ggplot(sres[sres$decay%in%c(100,100000)&sres$gamma%in%c(0.5,2),],aes(x=PC1,y=emissions,col=pop/totalPop))
g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',n=50)+xlab("PC1")+ylab("Normalized potential emissions")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_pc1-emissions_targeted.png'),width=30,height=25,units='cm')



g=ggplot(sres,aes(x=PC1,y=emissions,col=radius))
g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',span=0.5)+xlab("PC1")+ylab("Normalized potential emissions")+scale_color_continuous(name=expression(r[0]))+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_pc1-emissions_colr0.png'),width=30,height=25,units='cm')

g=ggplot(sres,aes(x=PC2,y=emissions,col=pop/totpop))
g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',n=50)+xlab("PC2")+ylab("Normalized potential emissions")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_pc2-emissions.png'),width=30,height=25,units='cm')

g=ggplot(sres,aes(x=PC2,y=emissions,col=radius))
g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',span=0.5)+xlab("PC2")+ylab("Normalized potential emissions")+scale_color_continuous(name=expression(r[0]))+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_pc2-emissions_colr0.png'),width=30,height=25,units='cm')

# same plots but relative to relative pop : as emissions are summed on clusters, does not make sense to look at absolute emissions
# RQ : before, not 'relative' but 'normalized' !
g=ggplot(sres,aes(x=PC1,y=emissions*totpop/pop,col=pop/totalPop))
g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',n=50)+xlab("PC1")+ylab("Relative potential emissions")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_pc1-relemissions.png'),width=30,height=25,units='cm')

g=ggplot(sres,aes(x=PC2,y=emissions*totpop/pop,col=pop/totalPop))
g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',n=50)+xlab("PC2")+ylab("Relative potential emissions")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_pc2-relemissions.png'),width=30,height=25,units='cm')


######


g=ggplot(sres,aes(x=PC1,y=efficiency,col=pop/totpop))
g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',n=50)+xlab("PC1")+ylab("1 - Normalized potential efficiency")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_pc1-efficiency.png'),width=30,height=25,units='cm')

g=ggplot(sres,aes(x=PC1,y=efficiency*totalPop/pop,col=pop/totpop))
g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',n=50)+xlab("PC1")+ylab("Relative potential unefficiency")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_pc1-relefficiency.png'),width=30,height=25,units='cm')


g=ggplot(sres,aes(x=PC2,y=efficiency,col=pop/totpop))
g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',n=50)+xlab("PC2")+ylab("1 - Normalized potential efficiency")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_pc2-efficiency.png'),width=30,height=25,units='cm')

g=ggplot(sres,aes(x=PC2,y=efficiency*totalPop/pop,col=pop/totpop))
g+geom_point()+facet_grid(decay~gamma)+geom_smooth(method='loess',n=50)+xlab("PC2")+ylab("Relative potential unefficiency")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_pc2-relefficiency.png'),width=30,height=25,units='cm')



g=ggplot(sres,aes(x=emissions,y=efficiency,col=PC1,size=pop/totalPop))
g+geom_point(alpha=0.5)+facet_grid(decay~gamma)+scale_color_viridis_c(name="PC1")+xlab("Normalized potential emissions")+ylab("1 - Normalized potential efficiency")+scale_size_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_emissions-efficiency_colpc1.png'),width=30,height=25,units='cm')

g=ggplot(sres,aes(x=emissions,y=efficiency,col=PC2,size=pop/totalPop))
g+geom_point(alpha=0.5)+facet_grid(decay~gamma)+scale_color_viridis_c(name="PC2")+xlab("Normalized potential emissions")+ylab("1 - Normalized potential efficiency")+scale_size_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_emissions-efficiency_colpc2.png'),width=30,height=25,units='cm')



##### pareto emissions / unefficiency
g=ggplot(sres,aes(x=1+emissions*totalPop/pop,y=1+efficiency*totalPop/pop,col=PC1,size=pop/totalPop))
g+geom_point(alpha=0.5)+scale_x_log10()+scale_y_log10()+facet_grid(decay~gamma)+scale_color_viridis_c(name="PC1")+xlab("Relative potential emissions")+ylab("Relative potential unefficiency")+scale_size_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_relemissions-relefficiency_colpc1_logscale.png'),width=30,height=25,units='cm')

## targeted
g=ggplot(sres[sres$gamma%in%c(0.5,2)&sres$decay%in%c(100,100000),],aes(x=1+emissions*totalPop/pop,y=1+efficiency*totalPop/pop,col=PC1,size=pop/totalPop))
g+geom_point(alpha=0.5)+scale_x_log10()+scale_y_log10()+facet_grid(decay~gamma)+
  scale_color_viridis_c(name="PC1")+xlab("Relative potential emissions")+ylab("Relative potential unefficiency")+
  scale_size_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_relemissions-relefficiency_colpc1_logscale_targeted.png'),width=30,height=25,units='cm')




g=ggplot(sres,aes(x=1+emissions*totalPop/pop,y=1+efficiency*totalPop/pop,col=PC2,size=pop/totalPop))
g+geom_point(alpha=0.5)+scale_x_log10()+scale_y_log10()+facet_grid(decay~gamma)+scale_color_viridis_c(name="PC2")+xlab("Relative potential emissions")+ylab("Relative potential unefficiency")+scale_size_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_morpho_relemissions-relefficiency_colpc2_logscale.png'),width=30,height=25,units='cm')




## total emissions as a function of morphology

g=ggplot(sres,aes(x=PC1,y=totalemissions,col=pop/totpop))
g+geom_point()+geom_smooth(method='loess',n=50)+xlab("PC1")+ylab("Total emissions")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_totemissions-pc1.png'),width=20,height=18,units='cm')

g=ggplot(sres,aes(x=PC2,y=totalemissions,col=pop/totpop))
g+geom_point()+xlab("PC2")+ylab("Total emissions")+scale_color_continuous(name="Relative\nPopulation")+stdtheme
ggsave(file=paste0(resdir,'aggreg_totemissions-pc2.png'),width=20,height=18,units='cm')


## try linear model to explain indicators as a function of morphology ?
# wont fit much given the ushaped relation
# : put square terms ?

fit = lm(emissions~moran+avgdist+entropy+slope,data=sres)
summary(fit)
AIC(fit)
# -> map the residuals ?

fitparams = lm(emissions~nwcol+popthq+nwthq+radius+gamma+decay,data=sres)
summary(fitparams)
AIC(fitparams) - AIC(fit)

summary(lm(emissions~paramrow,data=sres))


fitfull = lm(emissions~nwcol+popthq+nwthq+radius+gamma+decay+moran+avgdist+entropy+slope,data=sres)
summary(fitfull)
AIC(fitfull) - AIC(fit)
AIC(fitfull) - AIC(fitparams)
# -> fit improvement



# same with efficiency

fit = lm(efficiency~moran+avgdist+entropy+slope,data=sres)
summary(fit)
AIC(fit)
# -> map the residuals ?

fitparams = lm(efficiency~nwcol+popthq+nwthq+radius+gamma+decay,data=sres)
summary(fitparams)
AIC(fitparams) - AIC(fit)

summary(lm(efficiency~paramrow,data=sres))


fitfull = lm(efficiency~nwcol+popthq+nwthq+radius+gamma+decay+moran+avgdist+entropy+slope,data=sres)
summary(fitfull)
AIC(fitfull) - AIC(fit)
AIC(fitfull) - AIC(fitparams)

# rq : entropy switches from significant to non significant from emissions to efficiency

# add other indics ?
#

fitindics = lm(emissions~pophierarchy.rank+area+pop,data=sres)
summary(fitindics)

fitindics2 = lm(emissions~pophierarchy.rank+area+pop+avgpop+avgtotemissions,data=sres)
summary(fitindics2)
AIC(fitindics)-AIC(fitindics2)

g=ggplot(sres,aes(x=pop,y=emissions))
g+geom_point()+geom_smooth()
# -> explains the decreasing sign
g=ggplot(sres,aes(x=pop,y=efficiency))
g+geom_point()+geom_smooth()

fitindics = lm(efficiency~pophierarchy.rank+area+pop,data=sres)
summary(fitindics)

fitindics2 = lm(efficiency~pophierarchy.rank+area+pop+avgpop+avgtotemissions,data=sres)
summary(fitindics2)
AIC(fitindics)-AIC(fitindics2)


