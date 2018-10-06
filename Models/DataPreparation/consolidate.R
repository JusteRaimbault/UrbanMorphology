
##
# consolidate datasets


setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology'))

library(dplyr)
library(ggplot2)
source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

source('Models/StaticCorrelations/mapFunctions.R')

# network and pop indicators
areasize=100;offset=50;factor=0.5
indics = as.tbl(loadIndicatorData(paste0(Sys.getenv('CS_HOME'),"/UrbanMorphology/Data/StaticCorrelations/res/europecoupled_areasize",areasize,"_offset",offset,"_factor",factor,"_temp.RData")))


# different emissions gases
emissionsCO2 <- as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_C_2012/v432_CO2_excl_short-cycle_org_C_2012.txt',sep=';'))
emissionsCH4 <- as.tbl(read.csv('Data/EDGAR/v432_CH4_2012/v432_CH4_2012.txt',sep=';'))
emissionsN20 <- as.tbl(read.csv('Data/EDGAR/v432_N2O_2012/v432_N2O_2012.txt',sep=';'))
emissionsCO2short <- as.tbl(read.csv('Data/EDGAR/v432_CO2_org_short-cycle_C_2012/v432_CO2_org_short-cycle_C_2012.txt',sep=';'))

# summaries
summary(emissionsCO2$emission)
summary(emissionsCH4$emission)
summary(emissionsN20$emission)
summary(emissionsCO2short$emission)
# Lashof, D. A., & Ahuja, D. R. (1990). Relative contributions of greenhouse gas emissions to global warming. Nature, 344(6266), 529.
# => keep only CO2
# (for more thorough study should be more precise on ghg potential ?)

aggregems=data.frame(emissions=c(sum(emissionsCO2$emission),sum(emissionsCH4$emission),sum(emissionsN20$emission),sum(emissionsCO2short$emission)),
                    Type=c("CO2","CH4","N2O","CO2 (short)")
                     )
g <- ggplot(aggregems, aes(x="",y=emissions, fill=Type))
g+geom_bar(width = 1, stat = "identity")+xlab("")+ylab("Emissions (tons in 2012)")+stdtheme #+ coord_polar("y")
ggsave(file='Results/Stylized/allgases.png',width=18,height=15,units='cm')


cO2joined = left_join(emissionsCO2,emissionsCO2short,by=c('lat','lon'))
cO2joined[is.na(cO2joined)]=0
cO2 = data.frame(lat=cO2joined$lat,lon=cO2joined$lon,emission=cO2joined$emission.x+cO2joined$emission.y)
# aggregate by indic cells (requires data.R to have been run)
load('Models/DataPreparation/rowsemissions.RData')
cO2$indicrow = rows
cO2indics = cO2 %>% group_by(indicrow) %>% summarize(emission=sum(emission))
emission=cO2indics$emission;names(emission)=cO2indics$indicrow
# assumes indics has been loaded
indics$CO2 = emission[1:nrow(indics)]


# transportation emissions : short cycle for transport
em13aCDS <- as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_C_2012_IPCC_1A3a_CDS/v432_CO2_excl_short-cycle_org_C_2012_IPCC_1A3a_CDS.txt',sep=';'))
em13aCRS <- as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_C_2012_IPCC_1A3a_CRS/v432_CO2_excl_short-cycle_org_C_2012_IPCC_1A3a_CRS.txt',sep=';'))
em13aLTO <- as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_C_2012_IPCC_1A3a_LTO_IPCC_1A3a_LTO/v432_CO2_excl_short-cycle_org_C_2012_IPCC_1A3a_LTO_IPCC_1A3a_LTO.txt',sep=';'))
em13b <- as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_C_2012_IPCC_1A3b/v432_CO2_excl_short-cycle_org_C_2012_IPCC_1A3b.txt',sep=';'))
em13bshort <- as.tbl(read.csv('Data/EDGAR/v432_CO2_org_short-cycle_C_2012_IPCC_1A3b/v432_CO2_org_short-cycle_C_2012_IPCC_1A3b.txt',sep=';'))
em13ce <- as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_C_2012_IPCC_1A3c_1A3e/v432_CO2_excl_short-cycle_org_C_2012_IPCC_1A3c_1A3e.txt',sep=';'))
em13ceshort <- as.tbl(read.csv('Data/EDGAR/v432_CO2_org_short-cycle_C_2012_IPCC_1A3c_1A3e/v432_CO2_org_short-cycle_C_2012_IPCC_1A3c_1A3e.txt',sep=';'))
em13d <-  as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_2012_IPCC_1A3d_1C2/v432_CO2_excl_short-cycle_org_2012_IPCC_1A3d_1C2.txt',sep=';'))
em13dshort <- as.tbl(read.csv('Data/EDGAR/v432_CO2_org_short-cycle_C_2012_IPCC_1A3d_1C2/v432_CO2_org_short-cycle_C_2012_IPCC_1A3d_1C2.txt',sep=';'))

summary(em13aCDS$emission)
summary(em13aCRS$emission)

aggregems = data.frame(emissions = c(sum(em13aCDS$emission),sum(em13aCRS$emission),sum(em13aLTO$emission),sum(em13b$emission),sum(em13bshort$emission),
      sum(em13ce$emission),sum(em13ceshort$emission),sum(em13d$emission),sum(em13dshort$emission)),
    Type=c("Aviation (CD)","Aviation (CR)","Aviation (LT)","Road","Road (short)","Other","Other (short)","Shipping","Shipping (short)")
    )
g <- ggplot(aggregems, aes(x="",y=emissions, fill=Type))
g+geom_bar(width = 1, stat = "identity")+xlab("")+ylab("Emissions (tons in 2012)")+stdtheme #+ coord_polar("y")
ggsave(file='Results/Stylized/transport_C02.png',width=18,height=15,units='cm')


##
# same coordinates ? -> can merge with emissionsCO2, the most points for all CO2 datasets

transport=emissionsCO2
transport$emission=rep(0,nrow(transport))
for(type in list(em13aCDS,em13aCRS,em13aLTO,em13b,em13bshort,em13ce,em13ceshort,em13d,em13dshort)){
  cO2joined = left_join(transport,type,by=c('lat','lon'))
  cO2joined[is.na(cO2joined)]=0
  transport$emission=transport$emission+cO2joined$emission.y
}

# rq : share of transport ? 
# sum(transport$emission)/sum(cO2$emission)
# => 17%

transport$indicrow = rows
transportindics = transport %>% group_by(indicrow) %>% summarize(emission=sum(emission))
emission=transportindics$emission;names(emission)=transportindics$indicrow
# assumes indics has been loaded
indics$transport = emission[1:nrow(indics)]

# save consolidated indics
save(indics,file='Data/consolidated/indics.RData')






