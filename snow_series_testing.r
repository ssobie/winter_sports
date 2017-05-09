##Script to calculate the climdex indices from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file

##Updated version from compute.climdex.bccaq.r
##This computes all the climdex variables

source('/storage/data/projects/rci/assessments/code/snow.model.r',chdir=T)

library(ncdf4)
library(PCICt)

##---------------------------------------------------------------

snow.for.model <- function(var.names,gcm,scenario,interval,coords,
                           past.int,data.dir,write.dir) {

  pr.files <- list.files(path=paste(data.dir,gcm,sep=''),pattern='pr_gcm_prism',full.name=TRUE)
  pr.past.file <- pr.files[grep(past.int,pr.files)]
  run <- 'r1' ##strsplit(pr.past.file,'_')[[1]][11] ##[9]
  
  tasmax.files <- list.files(path=paste(data.dir,gcm,sep=''),pattern='tasmax_gcm_prism',full.name=TRUE)
  tasmax.past.file <- tasmax.files[grep(past.int,tasmax.files)]
  
  tasmin.files <- list.files(path=paste(data.dir,gcm,sep=''),pattern='tasmin_gcm_prism',full.name=TRUE)
  tasmin.past.file <- tasmin.files[grep(past.int,tasmin.files)]
  
  hist.dir <- paste(write.dir,gcm,'/',sep='')

  clim.all.files <- list.files(path=hist.dir,pattern='BCCAQ-PRISM',full.name=TRUE)
  clim.files <- clim.all.files[grep(scenario,clim.all.files)]
  clim.ncs <- lapply(clim.files,nc_open,write=TRUE)

  ##--------------------------------------------------------------
  print('Reading past')
  pr.past.nc <- nc_open(pr.past.file,write=FALSE)
  tasmax.past.nc <- nc_open(tasmax.past.file,write=FALSE)
  tasmin.past.nc <- nc_open(tasmin.past.file,write=FALSE)

  pr.scale <- 1
  temp.offset <- 0
  pr.units <- ncatt_get(pr.past.nc,'pr')$units
  if (pr.units == 'kg m-2 s-1')
    pr.scale <- 86400
  tx.units <- ncatt_get(tasmax.past.nc,'tasmax')$units
  if (tx.units == 'K')
    temp.scale <- 273
  print('Reading attributes')
  ##Attributes to retain
  lon <- ncvar_get(pr.past.nc,'lon')
  lat <- ncvar_get(pr.past.nc,'lat')  
  n.lon <- length(lon)
  n.lat <- length(lat)

  ##Combine the dates
  time.atts <- ncatt_get(pr.past.nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  
  pr.past.values <- ncvar_get(pr.past.nc,'time')
  pr.time.values <- pr.past.values
  pr.dates <- past.origin + pr.past.values*86400/24

  tasmax.past.values <- ncvar_get(tasmax.past.nc,'time')
  tasmax.time.values <- tasmax.past.values
  tasmax.dates <- as.Date(as.character(past.origin + tasmax.past.values*86400/24))

  tasmin.past.values <- ncvar_get(tasmin.past.nc,'time')
  tasmin.time.values <- tasmin.past.values

  ##--------------------------------------------------------------
  ##Compute climdex values and load into newly created climdex netcdf

  pr.past.subset <- ncvar_get(pr.past.nc,'pr',start=c(coords,1),count=c(1,1,-1))*pr.scale
  tasmax.past.subset <- ncvar_get(tasmax.past.nc,'tasmax',start=c(coords,1),count=c(1,1,-1))-temp.offset
  tasmin.past.subset <- ncvar_get(tasmin.past.nc,'tasmin',start=c(coords,1),count=c(1,1,-1))-temp.offset
      
  pr.subset <- pr.past.subset
  tasmax.subset <- tasmax.past.subset
  tasmin.subset <- tasmin.past.subset

  snow.objects <- snow.melt(pr.subset,tasmax.subset,tasmin.subset,tasmax.dates,lat[coords[2]])
  nc_close(pr.past.nc)
  nc_close(tasmax.past.nc)
  nc_close(tasmin.past.nc)
  return(snow.objects)

}

##**************************************************************************************

scenario <- 'rcp85'
past.int <- '1979-2016'
  
data.dir <-  '/storage/data/projects/rci/data/winter_sports/BCCAQ2/' 
write.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/' 
gcm <- 'ERA'

low.ix <- c(130,65)
high.ix <- c(141,106)
coords <- low.ix  
snow.objects <- snow.for.model(var.names=c('snowdepth','swe'),gcm,scenario,interval,coords,
                               past.int,data.dir,write.dir)




