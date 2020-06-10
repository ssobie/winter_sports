##Script to calculate the ARI (antecedent rainfall index)

ptm <- proc.time()

library(PCICt)
library(ncdf4)

source('/storage/home/ssobie/code/repos/winter_sports/format_gcm_prism_file_metadata.r')

##------------------------------------------------------------
##Pull the time series
time_series_from_nc <- function(nc) {

  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units

  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                     cal=time.calendar)
  time.values <- ncvar_get(nc,'time')
  time.series <- origin + time.values*86400

  time.min <- as.character(format(head(time.series,1),'%Y%m%d'))
  time.max <- as.character(format(tail(time.series,1),'%Y%m%d'))

  rv <- list(units=time.units,
             values=time.values,
             series=time.series,
             calendar=time.calendar,
             tmax=time.max,
             tmin=time.min)
  return(rv)
}

##------------------------------------------------------------

get_var_units <- function(var.name) {
   rv <- switch(var.name,
                pr='kg m-2 d-1',
                tasmax='degC',
                tasmin='degC')
   return(rv)
}

##------------------------------------------------------------
##Function From James
## n is the number of bytes (16 for short integer)

compute_scale_and_offset <- function(minimum, maximum, n) {
    # stretch/compress data to the available packed range
    scale.factor <- (maximum - minimum) / (2 ** n - 1)
    # translate the range to be symmetric about zero
    add.offset <- minimum + 2 ** (n - 1) * scale.factor
    c(scale.factor, add.offset)
}


##------------------------------------------------------------

get_scaling <- function(var.name) {
  ##Offset and scale for 16-bit integer
  if (grepl('tas',var.name)) { ##-100 to 100 degrees
    scale_factor <- 0.003051804
    add_offset <- 0.001525902
  }
  if (var.name=='pr') { ##0 to 1500 mm
    scale_factor <- 0.02288853
    add_offset <- 750.01144427
  }
  rv <- list(scale=scale_factor,
             offset=add_offset)
  return(rv)
}

##------------------------------------------------------------

pack_data_for_insertion <- function(input,var.name) {

  if (var.name == 'pr') {
     input[input >= 1500] <- 1490
  }
  scaling <- get_scaling(var.name)
  rv <- round((input - scaling$offset) / scaling$scale)
  rm(input)
  return(rv)

}


##------------------------------------------------------------

make_gcm_prism_file <- function(var.name,gcm,run,scenario,
                                gcm.file,ds.file,                                   
                                tmp.dir) {

  ds.nc <- nc_open(paste0(tmp.dir,ds.file))
  lon <- ncvar_get(ds.nc,'lon')
  nlon <- length(lon)
  lat <- ncvar_get(ds.nc,'lat')
  nlat <- length(lat)

  gcm.nc <- nc_open(paste0(tmp.dir,gcm.file))
  grid <- ncatt_get(gcm.nc,0)$grid_label
  ##Consider acquiring the other labels (gcm,scenario) similarly?
  time <- time_series_from_nc(gcm.nc)
  
  gcm.prism.file <- paste0(varname,'_day_BCCAQv2+PNWNAmet+PRISM_Vancouver_Whistler_',
                           gcm,'_historical+',scenario,'_',run,'_',
                           time$tmin,'-',time$tmax,'.nc')

  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat) 
  t.geog <- ncdim_def('time', time$units, time$values,
                       unlim=FALSE, calendar=time$calendar)

  var.geog <- ncvar_def(var.name, units=get_var_units(var.name),
                        dim=list(x.geog, y.geog, t.geog),
                        missval=32767.,prec='short')

  print('Creating netcdf')
  write.nc <- nc_create(paste0(tmp.dir,gcm.prism.file), var.geog,h_minfree=104857)

  print('Adding metadata')
  add_the_metadata(var.name,scenario,run,gcm.nc,write.nc)

  print('Adding Lon and Lat')
  ncvar_put(write.nc,'lon',lon)
  ncvar_put(write.nc,'lat',lat)
  nc_close(gcm.nc)

  rv <- list(file=gcm.prism.file,nc=write.nc,ds.nc=ds.nc)
  return(rv)
  
}

##------------------------------------------------------------

insert_file_data <- function(var.name,bccaq2.file,tmp.dir) {
 
  write.nc <- bccaq2.file$nc
  ds.nc <- bccaq2.file$ds.nc
  lon <- ncvar_get(ds.nc,'lon')
  nlon <- length(lon)
  lat <- ncvar_get(ds.nc,'lat')
  nlat <- length(lat)
  time <- ncvar_get(ds.nc,'time')
  ntime <- length(time)

  for (j in 1:nlat) {      
     print(paste0('Latitude: ',j,' of ',nlat))
     split.sub <- ncvar_get(ds.nc,var.name,start=c(1,j,1),count=c(-1,1,-1))
     split.insert <- pack_data_for_insertion(split.sub,var.name)

     ncvar_put(write.nc,var.name,split.insert,start=c(1,j,1),count=c(nlon,1,ntime))
     rm(split.sub)
   }
}

##------------------------------------------------------------

##************************************************************
testing <- FALSE


if (testing) {
   varname <- 'tasmax'
   gcm <- 'ACCESS1-0'
   run <- 'r1i1p1'
   tmpdir <- '/local_temp/ssobie/vw_meta/'
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
   }
}

scenario <- 'rcp85'
var.name <- varname
gcm.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/CMIP5/',gcm,'/')
gcm.files <- list.files(path=gcm.dir,pattern=scenario)
gcm.file <- gcm.files[grep(var.name,gcm.files)]
print(gcm.file)

bccaq2.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/',gcm,'/')
write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/vw_whistler_gcm_prism/',gcm,'/')
ds.file <- list.files(path=bccaq2.dir,pattern=paste0(var.name,'_gcm_prism_'))
print(ds.file)
if (length(ds.file) != 1) {
   stop('Incorrect number of downscaled files')
}

if (!file.exists(write.dir)) {
  dir.create(write.dir,recursive=T)
}

tmp.dir <- paste0(tmpdir,'combine_',gcm,'_',varname,'_',run,'_',scenario,'_tmp/')
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=T)
}
 
file.copy(from=paste0(gcm.dir,gcm.file),to=tmp.dir)
print('Done copying gcm file')
file.copy(from=paste0(bccaq2.dir,ds.file),to=tmp.dir)
print('Done copying DS file')
print('Done Copying')

##-----------------------------------------------

print('Making File') 

bccaqv2.file <- make_gcm_prism_file(var.name=varname,gcm=gcm,run=run,scenario=scenario,
                                    gcm.file=gcm.file,ds.file=ds.file,                              
                                    tmp.dir=tmp.dir)
file.remove(paste0(tmp.dir,gcm.file))

print('File data insertion') 
insert_file_data(var.name,bccaqv2.file,tmp.dir)
nc_close(bccaqv2.file$nc)
nc_close(bccaqv2.file$ds.nc)

print('Copying to write dir')
file.copy(from=paste0(tmp.dir,bccaqv2.file$file),to=write.dir,overwrite=TRUE)

file.remove(paste0(tmp.dir,bccaqv2.file$file))
rm.files <- list.files(path=paste0(tmp.dir,gcm,'/'),recursive=T,full.name=T)
file.remove(rm.files)
print('Copied back and cleaned up')

print('Elapsed time')
print(proc.time()-ptm)

