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
                   references="Dee, D.P. et al, (2011) The ERA‐Interim reanalysis: configuration and performance of the data assimilation system, Quarterly Journal of the Royal Meteorological Society, 137, 553-597")
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
##-------------------------------------------------------------------------------------

make_era_wind_file <- function(nc,new.dates,new.file,var.name) {

  time.units <- 'days since 1950-01-01'
  past.dates <- as.Date(as.character(format(netcdf.calendar(nc),'%Y-%m-%d')))
  time.calendar <- 'proleptic_gregorian'
  all.dates <- as.Date(format(new.dates,'%Y-%m-%d'))  #####c(past.dates,as.Date(format(new.dates,'%Y-%m-%d'))) 
  time.vals <- all.dates - as.Date('1950-01-01')

  ##--------------------------------------------------------------
  lon <- ncvar_get(nc,'longitude')
  lat <- ncvar_get(nc,'latitude')
  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, as.numeric(time.vals),
                        unlim=TRUE, calendar=time.calendar)

  var.geog <- ncvar_def(var.name, units='m s-1', dim=list(x.geog, y.geog, t.geog),
                        missval=-32768)
  file.nc <- nc_create(new.file, var.geog)

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

  ncvar_put(file.nc,varid='lon',vals=lon)
  ncvar_put(file.nc,varid='lat',vals=lat)

  ##Clear extraneous history
  ncatt_put(file.nc,varid=0,attname='history',attval='')

  return(file.nc)  
}

tmp.dir <- '/local_temp/ssobie/erai/'

if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=T)
}

data.dir <- paste0('/storage/data/projects/rci/data/winter_sports/ERA_INTERIM/wind/')
download.dir <- paste0(data.dir,'download/')

template.file <- 'uwind_12_2018-10-09.nc' ##Chosen at random for spatial info

file.copy(from=download.dir,to=tmp.dir,recursive=TRUE)
past.tmp.file <- paste0(tmp.dir,'download/',template.file)

nc <- nc_open(past.tmp.file)
past.lon <- ncvar_get(nc,'longitude')
past.lat <- ncvar_get(nc,'latitude')

new.dates <- seq(from=as.Date('1979-01-01'),by='day',to=as.Date('2018-12-31'))


uwind.file <- 'uwind_day_ERA-INTERIM_global_19790101-20181231.nc'
vwind.file <- 'vwind_day_ERA-INTERIM_global_19790101-20181231.nc'
wspd.file <- 'wspd_day_ERA-INTERIM_global_19790101-20181231.nc'

uwind.nc <- make_era_wind_file(nc,new.dates,paste0(tmp.dir,uwind.file),'uwind')
vwind.nc <- make_era_wind_file(nc,new.dates,paste0(tmp.dir,vwind.file),'vwind')
wspd.nc <- make_era_wind_file(nc,new.dates,paste0(tmp.dir,wspd.file),'wspd')

nlon <- nc$dim$longitude$len
nlat <- nc$dim$latitude$len
ntime <- length(new.dates)

nc_close(nc)

for (i in seq_along(new.dates)) {
  print(i)
  print(new.dates[i])
  ##------------------------------------
  ##Uwind
  uwind.zero.file <- paste0(tmp.dir,'download/uwind_00_',as.character(new.dates[i]),'.nc')
  zu.nc <- nc_open(uwind.zero.file)
  zuw <- ncvar_get(zu.nc,'u10')
  nc_close(zu.nc)

  uwind.twelve.file <- paste0(tmp.dir,'download/uwind_12_',as.character(new.dates[i]),'.nc')
  tu.nc <- nc_open(uwind.twelve.file)
  tuw <- ncvar_get(tu.nc,'u10')
  nc_close(tu.nc)
  uwind.mean <- apply(abind(zuw,tuw,along=3),c(1,2),mean)
  rm(tuw)
  rm(zuw)
  ncvar_put(uwind.nc,'uwind',uwind.mean,start=c(1,1,i),count=c(nlon,nlat,1))


  ##------------------------------------
  ##Vwind
  vwind.zero.file <- paste0(tmp.dir,'download/vwind_00_',as.character(new.dates[i]),'.nc')
  zv.nc <- nc_open(vwind.zero.file)
  zvw <- ncvar_get(zv.nc,'v10')
  nc_close(zv.nc)

  vwind.twelve.file <- paste0(tmp.dir,'download/vwind_12_',as.character(new.dates[i]),'.nc')
  tv.nc <- nc_open(vwind.twelve.file)
  tvw <- ncvar_get(tv.nc,'v10')
  nc_close(tv.nc)
  vwind.mean <- apply(abind(zvw,tvw,along=3),c(1,2),mean) 

  rm(tvw)
  rm(zvw)
  ncvar_put(vwind.nc,'vwind',vwind.mean,start=c(1,1,i),count=c(nlon,nlat,1))

  wspd.mean <- sqrt( (uwind.mean^2) + (vwind.mean^2))

  rm(uwind.mean)
  rm(vwind.mean)
  ncvar_put(wspd.nc,'wspd',wspd.mean,start=c(1,1,i),count=c(nlon,nlat,1))
  rm(wspd.mean)

}


##ncvar_put(vwind.nc,'uwind',vwind.mean)
##ncvar_put(wspd.nc,'wspd',wspd.mean)

nc_close(uwind.nc)
nc_close(vwind.nc)
nc_close(wspd.nc)

file.copy(from=paste0(tmp.dir,uwind.file),to=data.dir,overwrite=TRUE)
file.copy(from=paste0(tmp.dir,vwind.file),to=data.dir,overwrite=TRUE)
file.copy(from=paste0(tmp.dir,wspd.file),to=data.dir,overwrite=TRUE)

file.remove(paste0(tmp.dir,uwind.file))
file.remove(paste0(tmp.dir,vwind.file))
file.remove(paste0(tmp.dir,wspd.file))
file.remove(paste0(tmp.dir,'download/*nc'))