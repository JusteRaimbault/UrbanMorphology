
setwd(paste0(Sys.getenv('CS_HOME'),'/UrbanMorphology'))

library(dplyr)
library(ggplot2)
library(viridis)

#emissions <- as.tbl(read.csv('Data/EDGAR/v432_CH4_2012/v432_CH4_2012.txt',sep=';'))
emissions <- as.tbl(read.csv('Data/EDGAR/v432_CO2_excl_short-cycle_org_C_2012/v432_CO2_excl_short-cycle_org_C_2012.txt',sep=';'))

#substance = 'CH4'
substance = 'CO2'

year='2012'

emissions$log_emission = log(emissions$emission);emissions$log_emission[emissions$log_emission<0]=0

emissions$norm_log_emission = (emissions$log_emission - min(emissions$log_emission)) / (max(emissions$log_emission) - min(emissions$log_emission))

# test map with ggplot
g=ggplot(emissions,aes(x=lon,y=lat,fill=norm_log_emission))
g+geom_raster()+scale_fill_viridis(option = "magma", direction = -1)
ggsave(file=paste0('Results/Maps/Test/',substance,'_',year,'.png'),width=30,height=20,units='cm')




