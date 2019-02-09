##Script to merge the files for ERA Interim from 8 times daily to once daily

library(ncdf4)
library(PCICt)
library(abind)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R')

make_era_precip_file <- function(nc,new.dates,new.file) {

  time.units <- 'days since 1950-01-01'
  past.dates <- as.Date(as.character(format(netcdf.calendar(nc),'%Y-%m-%d')))
  time.calendar <- 'gregorian'
  all.dates <- c(past.dates,as.Date(format(new.dates,'%Y-%m-%d'))) 
  time.vals <- all.dates - as.Date('1950-01-01')

  ##--------------------------------------------------------------
  lon <- ncvar_get(nc,'longitude')
  lat <- ncvar_get(nc,'latitude')
  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, as.numeric(time.vals),
                        unlim=TRUE, calendar=time.calendar)

  var.geog <- ncvar_def('pr', units='m', dim=list(x.geog, y.geog, t.geog),
                        missval=-32768)
  file.nc <- nc_create(new.file, var.geog)
  ncvar_put(file.nc,varid='lon',vals=lon)
  ncvar_put(file.nc,varid='lat',vals=lat)
  return(file.nc)  
}

data.dir <- '/storage/data/projects/rci/data/winter_sports/ERA_INTERIM/pr/'
var.name <- 'pr'

##Existing Precipitation File
##Replace this with updated daily

past.file <- '/storage/data/projects/rci/data/winter_sports/ERA_INTERIM/pr/pr.nc'

nc <- nc_open(past.file)
past.lon <- ncvar_get(nc,'longitude')
past.lat <- ncvar_get(nc,'latitude')
new.dates <- seq(from=as.Date('2016-11-01'),by='day',to=as.Date('2018-10-31'))
new.file <- '/storage/data/projects/rci/data/winter_sports/ERA_INTERIM/pr/erai_pr_1979-2018.nc'

new.nc <- make_era_precip_file(nc,new.dates,new.file)

past.pr <- ncvar_get(nc,'pr')

new.pr <- array(0,c(nc$dim$longitude$len,nc$dim$latitude$len,length(new.dates)))

day.file <- paste0(data.dir,'daily/daily-tp_',as.character(new.dates[1]),'.nc')
dnc <- nc_open(day.file)
lon <- ncvar_get(dnc,'longitude')
lat <- ncvar_get(dnc,'latitude')
lon.ix <- which(lon %in% past.lon)
lon.st <- lon.ix[1]
lon.cnt <- length(lon.ix)
lat.ix <- which(lat %in% past.lat)
lat.st <- lat.ix[1]
lat.cnt <- length(lat.ix)
nc_close(dnc)

for (i in seq_along(new.dates)) {
  print(new.dates[i])
  day.file <- paste0(data.dir,'daily/daily-tp_',as.character(new.dates[i]),'.nc')
  dnc <- nc_open(day.file)
  dpr <- ncvar_get(dnc,'tp',start=c(lon.st,lat.st,1),count=c(lon.cnt,lat.cnt,-1))
  new.pr[,,i] <- dpr
  nc_close(dnc)
}

nc_close(nc)

all.pr <- abind(past.pr,new.pr)

ncvar_put(new.nc,'pr',all.pr)
nc_close(new.nc)