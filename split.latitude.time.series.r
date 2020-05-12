##Divide the BCCI file into latitude bands of 10 cells 

library(ncdf4)
library(PCICt)

ptm <- proc.time()

##----------------------------------------------------------------

time_series_from_nc <- function(nc) {

  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units

  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                     cal=time.calendar)
  time.values <- ncvar_get(nc,'time')
  time.series <- origin + time.values*86400

  rv <- list(units=time.units,
             values=time.values,
	     series=time.series,
             calendar=time.calendar)
  return(rv)
}

##-------------------------------------------------------------------

get_var_units <- function(var.name) {

  rv <- switch(var.name,
               pr='kg m-2 d-1',
	       tasmax='degC',
               tasmin='degC',
               swe='mm',
               snowdepth='m',
               scale='',
               slope='',
               freq='')
  return(rv)
}

##-----------------------------------------------------------------------

make_subset_file <- function(nc,var.name,lat.ix,split.dir,split.file) {

  bcci.time <- time_series_from_nc(nc)
  
  lon <- ncvar_get(nc,'lon')
  lat.sub <- ncvar_get(nc,'lat')[lat.ix]

  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat.sub)
  t.geog <- ncdim_def('time', bcci.time$units, bcci.time$values,
                       unlim=FALSE, calendar=bcci.time$calendar)

  var.geog <- ncvar_def(var.name, units=get_var_units(var.name), dim=list(x.geog, y.geog, t.geog),
                        missval=-32768)
  print('Making this file:')
  print(paste0(split.dir,split.file))
  file.nc <- nc_create(paste0(split.dir,split.file), var.geog) ##,h_minfree=104857)

  ncatt_put(file.nc,varid=var.name,attname='units',attval=get_var_units(var.name))
  ncatt_put(file.nc,varid=var.name,attname='_FillValue',attval=-32768)
  ncatt_put(file.nc,varid=var.name,attname='standard_name',attval=toupper(var.name))
  ncatt_put(file.nc,varid=var.name,attname='long_name',attval=toupper(var.name))

  ncatt_put(file.nc,varid='time',attname='units',attval=bcci.time$units)
  ncatt_put(file.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='calendar',attval=bcci.time$calendar)

  lon.atts <- ncatt_get(nc,'lon')
  lon.names <- names(lon.atts)
  for (j in 1:length(lon.atts))
    ncatt_put(file.nc,varid='lon',attname=lon.names[j],attval=lon.atts[[j]])

  lat.atts <- ncatt_get(nc,'lat')
  lat.names <- names(lat.atts)
  for (j in 1:length(lat.atts))
    ncatt_put(file.nc,varid='lat',attname=lat.names[j],attval=lat.atts[[j]])
  ncvar_put(file.nc,'lon',lon)
  ncvar_put(file.nc,'lat',lat.sub)

  nc_close(file.nc)
}


##------------------------------------------------------------------------------

split_into_lat_band <- function(bcci.nc,
		                  var.name,
		    		  lat.ix,lat.st,lat.en,lat.cnt,
		  		  split.file,
		    		  split.dir) {

  bcci.band <- ncvar_get(bcci.nc,var.name,start=c(1,lat.st,1),count=c(-1,lat.cnt,-1))
  make_subset_file(bcci.nc,var.name,lat.ix,split.dir,split.file)
  print('Adding data')
  sub.nc <- nc_open(paste0(split.dir,split.file),write=TRUE)
  ncvar_put(sub.nc,var.name,bcci.band)
  rm(bcci.band)
  nc_close(sub.nc)
  gc()
}

##------------------------------------------------------------

split_apart_vw_file <- function(var.name,gcm,input.file,
                                tmp.dir,write.dir) {

   nc <- nc_open(paste0(tmp.dir,input.file))
   split.dir <- paste0(tmp.dir,var.name,'_',gcm,'_split10/')
   if (!file.exists(split.dir)) {
      dir.create(split.dir,recursive=T)
   }

   ##Divide up the Van Whislter PRISM latitude (323 total) into one 13-cell and 31 10-cell bands:
   splits <- c(rep(1,13),rep(2:32,each=10))
   split.files <- rep('F',32)
   for (i in 1:32) {
      print(paste0('Band: ',i,' of 32'))      
      lat.ix <- which(splits %in% i)
      
      lat.st <- lat.ix[1]
      lat.en <- tail(lat.ix,1)
      lat.cnt <- diff(range(lat.ix))+1
  
      split.file <- gsub(paste0(var.name,'_'),
       	                 paste0(var.name,'_L',sprintf('%02d',i),'_',sprintf('%02d',lat.st),'-',sprintf('%02d',lat.en),'_'),
	     	         input.file)
      split.files[i] <- split.file
      split_into_lat_band(nc,var.name,
                          lat.ix,lat.st,lat.en,lat.cnt,
                          split.file,split.dir)
      ##print('Copying split file back from temp')
      ##print(split.file)
      ##file.copy(from=paste0(split.dir,split.file),to=write.dir,overwrite=TRUE)                   
   }
   print('Copy files to write dir')
   file.copy(from=split.dir,to=write.dir,recursive=TRUE,overwrite=TRUE)
   files.rm <- list.files(path=split.dir,full.name=T)
   file.remove(files.rm)
   nc_close(nc)    
   file.remove(paste0(tmp.dir,input.file))
   print('Done with cleanup')
   ###Modify to return the split directory and a list of the split files
   ##rv <- list(dir=.dir,files=split.files)
   ##return(rv)

}

gcm <- 'ERA5'
tmpdir <- '/local_temp/ssobie'

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}

model <- gcm
tmp.dir <- paste0(tmpdir,'/',gcm,'_snow_split/')

if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir)
}

##--------------------------------------------------------------------------------------
##Downscaled GCM simulations

if (1==0) { 
  ##Empty snow files
  interval <- '19500101-21001231'
  read.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/templates/',model,'/')
  write.dir <- read.dir
  swe.file <- paste0("swe_BCCAQ2-PRISM_",model,"_",interval,".nc")
  file.copy(from=paste0(read.dir,swe.file),to=tmp.dir,overwrite=TRUE)
  snow.file <- paste0("snowdepth_BCCAQ2-PRISM_",model,"_",interval,".nc")
  file.copy(from=paste0(read.dir,snow.file),to=tmp.dir,overwrite=TRUE)
  print('Done copying snow files')    
  split_apart_vw_file('swe',model,swe.file,tmp.dir,write.dir)
  split_apart_vw_file('snowdepth',model,snow.file,tmp.dir,write.dir)
  file.remove(paste0(tmp.dir,snow.file))
  file.remove(paste0(tmp.dir,swe.file))
}

if (1==0) { 
  ##Input GCM-PRISM files
  read.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/',model,'/')
  write.dir <- read.dir
  interval <- '1950-2100'
  pr.file <- paste0("pr_gcm_prism_",model,"_VW_",interval,".nc")
  file.copy(from=paste0(read.dir,pr.file),to=tmp.dir,overwrite=TRUE)
  split_apart_vw_file('pr',model,pr.file,tmp.dir,write.dir)

  tasmax.file <- paste0("tasmax_gcm_prism_",model,"_VW_",interval,".nc")
  file.copy(from=paste0(read.dir,tasmax.file),to=tmp.dir,overwrite=TRUE)
  split_apart_vw_file('tasmax',model,tasmax.file,tmp.dir,write.dir)

  tasmin.file <- paste0("tasmin_gcm_prism_",model,"_VW_",interval,".nc")
  file.copy(from=paste0(read.dir,tasmin.file),to=tmp.dir,overwrite=TRUE)
  split_apart_vw_file('tasmin',model,tasmin.file,tmp.dir,write.dir)

  file.remove(paste0(tmp.dir,pr.file))
  file.remove(paste0(tmp.dir,tasmax.file))
  file.remove(paste0(tmp.dir,tasmin.file))

}


##--------------------------------------------------------------------------------------
##Reanalysis and Gridded observations

##Calibration Parameter files
if (1==1) {
  read.dir <- '/storage/data/projects/rci/data/winter_sports/'
  write.dir <- read.dir
  scale.file <- paste0("scale_hyper_snow_calibrated_parameter_",model,"_prism_TPS.nc")
  file.copy(from=paste0(read.dir,scale.file),to=tmp.dir,overwrite=TRUE)
  slope.file <- paste0("slope_hyper_snow_calibrated_parameter_",model,"_prism_TPS.nc")
  file.copy(from=paste0(read.dir,slope.file),to=tmp.dir,overwrite=TRUE)
  freq.file <- paste0("freq_hyper_snow_calibrated_parameter_",model,"_prism_TPS.nc")
  file.copy(from=paste0(read.dir,freq.file),to=tmp.dir,overwrite=TRUE)

  split_apart_vw_file('scale',model,scale.file,tmp.dir,write.dir)
  split_apart_vw_file('slope',model,slope.file,tmp.dir,write.dir)
  split_apart_vw_file('freq',model,freq.file,tmp.dir,write.dir)
}


  ##Empty Snow files

if (1==0) { 
  interval <- '1945-2012'
  read.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/',model,'/')
  write.dir <- read.dir
  swe.file <- paste0("swe_BCCAQ2-PRISM_",model,"_",interval,".nc")
  file.copy(from=paste0(read.dir,swe.file),to=tmp.dir,overwrite=TRUE)
  snow.file <- paste0("snowdepth_BCCAQ2-PRISM_",model,"_",interval,".nc")
  file.copy(from=paste0(read.dir,snow.file),to=tmp.dir,overwrite=TRUE)
  print('Done copying snow files')    
  split_apart_vw_file('swe',model,swe.file,tmp.dir,write.dir)
  split_apart_vw_file('snowdepth',model,snow.file,tmp.dir,write.dir)
  file.remove(paste0(tmp.dir,snow.file))
  file.remove(paste0(tmp.dir,swe.file))
}

if (1==0) { 
  ##Input GCM-PRISM files
  read.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/',model,'/')
  write.dir <- read.dir
  interval <- '19450101-20121231'
  pr.file <- paste0("pr_gcm_prism_",model,"_VW_",interval,".nc")
  file.copy(from=paste0(read.dir,pr.file),to=tmp.dir,overwrite=TRUE)
  split_apart_vw_file('pr',model,pr.file,tmp.dir,write.dir)

  tasmax.file <- paste0("tasmax_gcm_prism_",model,"_VW_",interval,".nc")
  #file.copy(from=paste0(read.dir,tasmax.file),to=tmp.dir,overwrite=TRUE)
  #split_apart_vw_file('tasmax',model,tasmax.file,tmp.dir,write.dir)

  tasmin.file <- paste0("tasmin_gcm_prism_",model,"_VW_",interval,".nc")
  #file.copy(from=paste0(read.dir,tasmin.file),to=tmp.dir,overwrite=TRUE)
  #split_apart_vw_file('tasmin',model,tasmin.file,tmp.dir,write.dir)

  file.remove(paste0(tmp.dir,pr.file))
#  file.remove(paste0(tmp.dir,tasmax.file))
#  file.remove(paste0(tmp.dir,tasmin.file))

}

