##Script to pull out the time series of pr, tasmax, tamsin from the driving models 
library(ncdf4)
library(rgdal)
library(raster)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##----------------------------------------------------------------------------
get_max_mask <- function(swe.file,swe.dir,clip.shp) {

   swe.data <- brick(paste0(swe.dir,swe.file))
   swe.max <- calc(swe.data,max,na.rm=T)
   swe.clip <- mask(swe.max,clip.shp)
   swe.mask <- swe.clip*0 + 1
   swe.mask[swe.clip > 5] <- NA
   
   return(swe.mask)
}



area_average_swe_series <- function(gcm,var.name,swe.file,swe.dir,clip.shp,swe.clip.mask) {

   swe.data <- brick(paste0(swe.dir,swe.file))

   nc <- nc_open(paste0(swe.dir,swe.file))
   var.time <- netcdf.calendar(nc)
   nc_close(nc)
   var.years <- format(var.time,'%Y')
   years <- unique(var.years)

   swe.series <- rep(0,length(var.years))
   for (y in seq_along(years)) { 
      print(years[y])
      time.ix <- grep(years[y],var.years)
      swe.sub <- subset(swe.data,time.ix)      
      swe.mask <- mask(swe.sub,clip.shp) * swe.clip.mask
      swe.series[time.ix] <- round(as.numeric(cellStats(swe.mask,mean,na.rm=T))*1000)
   }

  rv <- list(time=format(var.time,'%Y-%m-%d'),series=swe.series)

  return(rv)
}


##----------------------------------------------------------------------------
##---------------------------------------------------------------------
read.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/'
##gcm.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
##              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

##gcm <- 'ACCESS1-0'
##tmpdir <- '/local_temp/ssobie/reg/'

tmp.dir <- paste0(tmpdir,'/regional_',gcm,'/')
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=T)
}

var.name <- 'swe'
shp <- readOGR('/storage/data/projects/rci/data/assessments/metro_van/shapefiles/','MVWaterSheds', stringsAsFactors=F)
clip.shp <- spTransform(shp,CRS("+init=epsg:4326"))


##for (i in seq_along(gcm.list)) {
##   gcm <- gcm.list[i]
   print(gcm)
   swe.dir <- paste0(read.dir,'calibrated_',gcm,'_PNWNAmet_prism_tps/')

   ##Find the permanently snow cover cells and exclude them from the average
   swe.max.file <- list.files(path=swe.dir,pattern=paste0(var.name,'_annual_maximum_BCCAQ2-PRISM'))
   swe.max.mask <- get_max_mask(swe.max.file,swe.dir,clip.shp)

   swe.file <- list.files(path=swe.dir,pattern=paste0(var.name,'_BCCAQ2-PRISM'))

   file.copy(from=paste0(swe.dir,swe.file),to=tmp.dir,overwrite=TRUE)

   swe.info <-  area_average_swe_series(gcm,var.name,swe.file,tmp.dir,clip.shp,swe.max.mask)

   output <- cbind(swe.info$time,swe.info$series)
   output <- rbind(c('Dates','SWE'),output)
   write.table(output,
               file=paste('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/',
                           'mv_watersheds_',gcm,'_PNWNAmet_PRISM_TPS_snow_model_data_masked.csv',sep=''),
                           sep=',',row.name=FALSE,col.name=FALSE,quote=FALSE)

   file.remove(paste0(tmp.dir,swe.file))

##browser()
##  }



