##Script to convert the MODIS tiff files into netcdf for easier use

library(raster)
library(ncdf4)


##-------------------------------------------------------------------------------------
##Global ERA-Interim Attributes
get_global_atts <- function() {
  global.atts <- list(institution="National Snow and Ice Data Center",
                   contact="NSIDC",
                   Conventions="CF-1.6",
                   institute_id ="NSIDC",
                   domain='US / Southern Canada',
                   creation_date=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%Z'),
                   frequency='day',
                   product="reanalysis",
                   modeling_realm="surface",
                   project_id='SNODAS',
                   references="National Operational Hydrologic Remote Sensing Center. 2004. Snow Data Assimilation System (SNODAS) Data Products at NSIDC, Version 1. [SWE]. Boulder, Colorado USA. NSIDC: National Snow and Ice Data Center. doi: https://doi.org/10.7265/N5TB14TC. 12-Feb-2019.")
}


##-------------------------------------------------------------------------------------
get_standard_atts <- function(var.name) {
  lon.atts <- list(standard_name="longitude",long_name = "longitude",
                   units = "degrees_east",axis = "X")

  lat.atts <- list(standard_name="latitude",long_name = "latitude",
                   units = "degrees_north",axis = "Y")

  swe.atts <- list(standard_name = "snow_water_equivalent",
                  long_name = "Snow Water Equivalent",
                  missing_value = 1.e+20,
                  cell_methods = "time: mean")

  var.atts <- swe.atts ##switch(var.name,
  rv <- list(lon=lon.atts,
             lat=lat.atts,
             var=var.atts)
  return(rv)              
}
##-------------------------------------------------------------------------------------

##Function to create full time series file


make_snodas_netcdf_series <- function(var.name,day.dates,       
                                      base.file,write.file,
                                      lon.ix,lat.ix,
                                      tmp.dir) {
 
  bnc <- nc_open(paste0(tmp.dir,base.file))
  blon <- ncvar_get(bnc,'lon')
  blat <- ncvar_get(bnc,'lat')

  lon <- blon##[lon.ix]
  lat <- blat##[lat.ix]
 
  time.units <- 'days since 1950-01-01'
  time.vals <- day.dates - as.Date('1950-01-01')
  time.calendar <- 'gregorian'


  ##--------------------------------------------------------------
  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, as.numeric(time.vals),
                        unlim=FALSE, calendar=time.calendar)

  var.geog <- ncvar_def(var.name, units='mm', dim=list(x.geog, y.geog, t.geog),
                        missval=-32768)
  file.nc <- nc_create(paste0(tmp.dir,write.file), var.geog)

  ncatt_put(file.nc,varid='time',attname='units',attval=time.units)
  ncatt_put(file.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='calendar',attval=time.calendar)

  standard.atts <- get_standard_atts(var.name)
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

  print('Global atts')
  ##Global Attributes
  global.atts <- get_global_atts()
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(file.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  ncvar_put(file.nc,'lon',lon)
  ncvar_put(file.nc,'lat',lat)

  ##Clear extraneous history
  ncatt_put(file.nc,varid=0,attname='history',attval='')


  nc_close(bnc)
  nc_close(file.nc)
}

##-------------------------------------------------------------------------------------

add_data_to_netcdf <- function(var.name,lon.ix,lat.ix,tst,tct,
                               in.file,write.file,wnc,tmp.dir) {

  inc <- nc_open(paste0(tmp.dir,in.file))                              
  sub.data <- ncvar_get(inc,var.name)##[lon.ix,lat.ix,]
  ncvar_put(wnc,varid=var.name,vals=sub.data,start=c(1,1,tst),count=c(-1,-1,tct))
  nc_close(inc)
}

##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------

var.name <- 'swe'

##Subset for Van-Whistler
##lon.ix <- 717:1319
##lat.ix <- 2874:3350

lon.ix <- 1:8192
lat.ix <- 1:4096


proj.dir <- '/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/'
write.dir <- proj.dir
tmp.dir <- '/local_temp/ssobie/snodas/'

if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}

day.dates <- seq(from=as.Date(paste0('2010-01-01')),by='day',to=as.Date(paste0('2018-12-31')))

base.file <- paste0(var.name,'_snodas_unmasked_20100101-20100131.nc')
file.copy(from=paste0(proj.dir,base.file),to=tmp.dir,overwrite=T)

##write.file <- paste0(var.name,'_snodas_van_whistler_20100101-20181231.nc')
write.file <- paste0(var.name,'_snodas_us_canada_20100101-20181231.nc')

make_snodas_netcdf_series(var.name,day.dates,       
                          base.file,write.file,
                          lon.ix,lat.ix,
                          tmp.dir)

##browser()


month.dates <- seq(from=as.Date(paste0('2010-01-01')),by='month',to=as.Date(paste0('2018-12-31')))

ncdf.files <- list.files(path=proj.dir,pattern=var.name)

wnc <- nc_open(paste0(tmp.dir,write.file),write=TRUE)

for (mn in seq_along(month.dates)) {
  print(month.dates[mn])
 
  in.file <- ncdf.files[grep(format(month.dates[mn],'%Y%m'),ncdf.files)]

  tix <- grep(format(month.dates[mn],'%Y-%m'),day.dates)
  tst <- tix[1]
  tct <- length(tix)
  mn.day.dates <- day.dates[tix]

  file.copy(from=paste0(proj.dir,in.file),to=tmp.dir,overwrite=T)

  add_data_to_netcdf(var.name,lon.ix,lat.ix,tst,tct,
                     in.file,write.file,wnc,tmp.dir)
  file.remove(paste0(tmp.dir,in.file))
}

nc_close(wnc)

file.copy(from=paste0(tmp.dir,write.file),to=proj.dir,overwrite=T)
