##Script to merge the files for ERA Interim from 8 times daily to once daily

library(ncdf4)
library(PCICt)
library(abind)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R')


##-------------------------------------------------------------------------------------
##Global ERA-Interim Attributes
get_global_atts <- function() {
  global.atts <- list(institution="European Centre for Medium-Range Weather Forecasts",
                   contact="ECMWF",
                   Conventions="CF-1.6",
                   institute_id ="ECMWF",
                   domain='Global',
                   creation_date=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%Z'),
                   frequency='day',
                   product="reanalysis",
                   modeling_realm="atmos",
                   project_id='ERA-Interim',
                   references="Dee, D.P. et al, (2011) The ERA-Interim reanalysis: configuration and performance of the data assimilation system, Quarterly Journal of the Royal Meteorological Society, 137, 553-597")
}

##-------------------------------------------------------------------------------------
get_standard_atts <- function(var.name) {
  lon.atts <- list(standard_name="longitude",long_name = "longitude",
                   units = "degrees_east",axis = "X")

  lat.atts <- list(standard_name="latitude",long_name = "latitude",
                   units = "degrees_north",axis = "Y")

  pr.atts <- list(standard_name = "total_precipitation",
                  long_name = "Precipitation",
                  missing_value = 1.e+20,
                  cell_methods = "time: sum")
  tas.atts <- list(standard_name = "air_temperature",
                   missing_value = 1.e+20)
  uwind.atts <- list(standard_name = "10 metre U wind component",
                   missing_value = 1.e+20)
  vwind.atts <- list(standard_name = "10 metre U wind component",
                   missing_value = 1.e+20)

  var.atts <- switch(var.name,
                     pr=pr.atts,
                     tas=tas.atts,
                     sp=sp.atts,
                     rhs=rhs.atts,
                     huss=huss.atts,
                     uwind=uwind.atts,
                     vwind=vwind.atts)

  rv <- list(lon=lon.atts,
             lat=lat.atts,
             var=var.atts)
  return(rv)

}

##-------------------------------------------------------------------------------------
get_variable_units <- function(var.name) {

  units <- list(pr="kg m-2 day-1",rhs='%',huss='kg kg-1',
                tasmax='degC',tasmin='degC',tasday='degC',tashour='degC',tasrange='degC',
                tasskew='',sp='Pa',uwind='m s-1',vwind='m s-1',dewpoint='degC',
                tcc='fraction',lcc='fraction',rain='kg m-2')
  return(units[[var.name]])
}

##-------------------------------------------------------------------------------------
get_variable_specific_atts <- function(var.name) {

  pr.day.atts <- list(units = get_variable_units('pr'))
  tasmax.atts <- list(long_name = "Daily Maximum Near-Surface Air Temperature",
                      cell_methods = "time: maximum",units=get_variable_units('tasmax'))
  tasmin.atts <- list(long_name = "Daily Minimum Near-Surface Air Temperature",
                      cell_methods = "time: minimum",units=get_variable_units('tasmin'))
  tas.atts <- list(long_name = "Daily Average Near-Surface Air Temperature",
                      cell_methods = "time: mean",units=get_variable_units('tas'))
  uwind.atts <- list(long_name = "10 metre U wind component",
                     units=get_variable_units('uwind'))
  vwind.atts <- list(long_name = "10 metre V wind component",
                     units=get_variable_units('vwind'))

  var.atts <- switch(var.name,
                     pr=pr.day.atts,
                     tasmax=tasmax.atts,
                     tasmin=tasmin.atts,
                     tas=tas.atts,
                     uwind=uwind.atts,
                     vwind=vwind.atts)

  rv <- list(var=var.atts)
  return(rv)
}



make_era_tas_file <- function(nc,new.dates,new.file,var.name) {

  time.units <- 'days since 1950-01-01'
  past.dates <- as.Date(as.character(format(netcdf.calendar(nc),'%Y-%m-%d')))
  time.calendar <- 'gregorian'
  all.dates <- c(past.dates,as.Date(format(new.dates,'%Y-%m-%d'))) 
  time.vals <- all.dates - as.Date('1950-01-01')

  ##--------------------------------------------------------------
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, as.numeric(time.vals),
                        unlim=TRUE, calendar=time.calendar)

  var.geog <- ncvar_def(var.name, units='degC', dim=list(x.geog, y.geog, t.geog),
                        missval=-32768)
  file.nc <- nc_create(new.file, var.geog)

  ncvar_put(file.nc,varid='lon',vals=lon)
  ncvar_put(file.nc,varid='lat',vals=lat)

  ncatt_put(file.nc,varid='time',attname='units',attval=time.units)
  ncatt_put(file.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='calendar',attval=time.calendar)

  standard.atts <- get_standard_atts(var.name)
  variable.atts <- get_variable_specific_atts(var.name)
  print('Lon names')
  lon.names <- names(standard.atts$lon)
  for (j in 1:length(standard.atts$lon))
    ncatt_put(file.nc,varid='lon',attname=lon.names[j],attval=standard.atts$lon[[j]])
  print('Lat names')
  lat.names <- names(standard.atts$lat)
  for (j in 1:length(standard.atts$lat))
    ncatt_put(file.nc,varid='lat',attname=lat.names[j],attval=standard.atts$lat[[j]])

  print('Standard names')
  var.names <- names(standard.atts$var)
  for (j in 1:length(standard.atts$var))
    ncatt_put(file.nc,varid=var.name,attname=var.names[j],attval=standard.atts$var[[j]])
  print('Variable names')
  var.names <- names(variable.atts$var)
  for (j in 1:length(variable.atts$var))
    ncatt_put(file.nc,varid=var.name,attname=var.names[j],attval=variable.atts$var[[j]])

  print('Global atts')
  ##Global Attributes
  global.atts <- get_global_atts()
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(file.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  ##Clear extraneous history
  ncatt_put(file.nc,varid=0,attname='history',attval='')

  return(file.nc)  
}

past.file <- 'tasmin_day_ERA-INTERIM_global_19790101-20181031.nc'
var.name <- 'tasmin'
era.var <- 'mn2t'
tmp.dir <- '/local_temp/ssobie/erai/'

if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=T)
}

data.dir <- paste0('/storage/data/projects/rci/data/winter_sports/ERA_INTERIM/',var.name,'/')
download.dir <- paste0(data.dir,'download/')

##past.file <- paste0('/storage/data/projects/rci/data/winter_sports/ERA_INTERIM/',var.name,'/',var.name,'.nc')
past.dir <- paste0('/storage/data/projects/rci/data/winter_sports/ERA_INTERIM/',var.name,'/')


file.copy(from=paste0(past.dir,past.file),to=tmp.dir)
file.copy(from=download.dir,to=tmp.dir,recursive=TRUE)
past.tmp.file <- paste0(tmp.dir,past.file)

nc <- nc_open(past.tmp.file)
past.lon <- ncvar_get(nc,'lon')
past.lat <- ncvar_get(nc,'lat')
new.dates <- seq(from=as.Date('2018-11-01'),by='day',to=as.Date('2018-12-31'))

##new.file <- paste0('/storage/data/projects/rci/data/winter_sports/ERA_INTERIM/',var.name,'/erai_',var.name,'_1979-2018.nc')
new.file <- paste0('era_interim_global_',var.name,'_19790101-20181231.nc')

new.nc <- make_era_tas_file(nc,new.dates,paste0(tmp.dir,new.file),var.name)

past.tas <- ncvar_get(nc,var.name)

new.tas <- array(0,c(nc$dim$lon$len,nc$dim$lat$len,length(new.dates)))

day.file <- paste0(tmp.dir,'download/',var.name,'_00_',as.character(new.dates[1]),'.nc')
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
  zero.file <- paste0(tmp.dir,'download/',var.name,'_00_',as.character(new.dates[i]),'.nc')
  znc <- nc_open(zero.file)
  ztas <- ncvar_get(znc,era.var,start=c(lon.st,lat.st,1),count=c(lon.cnt,lat.cnt,-1))
  nc_close(znc)
  twelve.file <- paste0(tmp.dir,'download/',var.name,'_12_',as.character(new.dates[i]),'.nc')
  tnc <- nc_open(twelve.file)
  ttas <- ncvar_get(tnc,era.var,start=c(lon.st,lat.st,1),count=c(lon.cnt,lat.cnt,-1))
  nc_close(tnc)
  fx <- switch(var.name,tasmax=max,tasmin=min)

  tas <- apply(abind(ztas,ttas,along=3),c(1,2),fx)
  new.tas[,,i] <- tas
  rm(ttas)
  rm(ztas)
  rm(tas)

}
nc_close(nc)

all.tas <- abind(past.tas,new.tas,along=3)

ncvar_put(new.nc,var.name,all.tas)
nc_close(new.nc)

file.copy(from=paste0(tmp.dir,new.file),to=data.dir,overwrite=TRUE)