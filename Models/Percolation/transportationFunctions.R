
library(sf)
library(igraph)

source('percolationFunctions.R')

#'
#'
#'
transportEmissionsMSE <- function(indics,gamma,decay,xcol="lonmin",ycol="latmin",popcol="totalPop"){
  # generate flows
  sppoints <- st_as_sf(d, coords = c(xcol, ycol), crs = 4326) %>% st_transform(3035)
  flows = interactionPotential(sppoints,gamma,decay)

  # assignate : must first construct the graph : aggregate at this stage ?
  
  #graph_from_adj_list(,mode = 'all')
  
  
  
}
