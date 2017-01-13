##Script to convert the MODIS tiff files into netcdf for easier use

library(ncdf4)

##Function to create month file

aqua.dir <- '/storage/data/projects/rci/data/winter_sports/MODIS_AQUA_NETCDF/'
terra.dir <- '/storage/data/projects/rci/data/winter_sports/MODIS_TERRA_NETCDF/'


year.dates <- seq(from=as.Date('2002-01-01'),by='month',to=as.Date('2015-12-31'))

yrs <- format(year.dates,'%Y')
mns <- format(year.dates,'%m')

dates <- paste0(yrs,mns)

prefix <- 'snc.modis.terra.'

for (i in seq_along(dates)) {

  aqua.file <- list.files(path=aqua.dir,pattern=dates[i],full.name=TRUE)
  terra.file <- list.files(path=terra.dir,pattern=dates[i],full.name=TRUE)
  merged.file <- gsub(pattern='terra',replacement='merged',list.files(path=terra.dir,pattern=dates[i]))
  copy.terra <- paste0('cp ',terra.file,' /storage/data/projects/rci/data/winter_sports/MODIS_MERGED/',merged.file)
  system(copy.terra)

  aqua.nc <- nc_open(aqua.file)
  terra.nc <- nc_open(terra.file)

  aqua.data <- ncvar_get(aqua.nc,'snc')
  aqua.flag <- !is.na(aqua.data)
  terra.data <- ncvar_get(terra.nc,'snc')
  terra.flag <- !is.na(terra.data)
  nc_close(aqua.nc)
  nc_close(terra.nc)

  merged.data <- array(NA,dim(terra.data))
  merged.data[aqua.flag] <- aqua.data[aqua.flag]
  merged.data[terra.flag] <- terra.data[terra.flag]
  
  aqua.snow  <- which(aqua.data < 250) 
  merged.snow <- which(merged.data < 250) 
      
  if (length(aqua.snow) > 0) {
    merged.check <- !(aqua.snow %in% merged.snow)
    if (sum(merged.check) > 0) {
       merged.add <- aqua.snow[merged.check]
       merged.data[merged.add] <- aqua.data[merged.add]
    }
  }
      
  merged.nc <- nc_open(paste0('/storage/data/projects/rci/data/winter_sports/MODIS_MERGED/',merged.file),write=TRUE)

  ncvar_put(merged.nc,'snc',merged.data)
  nc_close(merged.nc)

}