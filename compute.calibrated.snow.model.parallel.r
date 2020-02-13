##Script to calculate the climdex indices from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file

##Updated version from compute.climdex.bccaq.r
##This computes all the climdex variables

source('/storage/home/ssobie/code/repos/winter_sports/snow.model.calibrated.r',chdir=T)

library(ncdf4)
library(PCICt)
library(zoo)

library(doParallel)
registerDoParallel(cores=4) # or some other number that you're comfortable with.


##---------------------------------------------------------------
separate_into_list <- function(data.subset) {
    return(lapply(seq_len(nrow(data.subset)), function(k) data.subset[k,]))
}


snow_for_model <- function(var.names,
                           pr.file,tasmax.file,tasmin.file,snow.file,swe.file,
                           data.dir,write.dir,tmp.dir) {

  atm <- proc.time()

  ##----------------------
  ##Slopes and Aspects from PRISM
  as.dir <- '/storage/data/projects/rci/data/prism/'
  slopes.nc <- nc_open(paste0(as.dir,'prism_slopes.nc'))
  bc.slopes <- ncvar_get(slopes.nc,'Band1')/90*pi/2
  bc.lon <- ncvar_get(slopes.nc,'lon')
  bc.lat <- ncvar_get(slopes.nc,'lat')
  nc_close(slopes.nc)

  aspects.nc <- nc_open(paste0(as.dir,'prism_aspects.nc'))
  bc.aspects <- ncvar_get(aspects.nc,'Band1')/360*2*pi
  nc_close(aspects.nc)
  ##----------------------
  ##Calibration Parameters
  cal.dir <- '/storage/data/projects/rci/data/winter_sports/'
  scale.nc <- nc_open(paste0(cal.dir,'scale_hyper_snow_calibrated_parameter.nc'))
  cal.scales <- ncvar_get(scale.nc,'scale')
  cal.lon <- ncvar_get(scale.nc,'lon')
  cal.lat <- ncvar_get(scale.nc,'lat')
  nc_close(scale.nc)

  slope.nc <- nc_open(paste0(cal.dir,'slope_hyper_snow_calibrated_parameter.nc'))
  cal.slopes <- ncvar_get(slope.nc,'slope')
  nc_close(slope.nc)

  freq.nc <- nc_open(paste0(cal.dir,'freq_hyper_snow_calibrated_parameter.nc'))
  cal.freqs <- ncvar_get(freq.nc,'freq')
  nc_close(freq.nc)

  ##----------------------
  ##Blank snow files
  depth.nc <- nc_open(paste0(tmp.dir,snow.file),write=TRUE)
  swe.nc <- nc_open(paste0(tmp.dir,swe.file),write=TRUE)

  ##--------------------------------------------------------------
  print('Reading past')
  pr.nc <- nc_open(paste0(tmp.dir,pr.file),write=FALSE)
  tasmax.nc <- nc_open(paste0(tmp.dir,tasmax.file),write=FALSE)
  tasmin.nc <- nc_open(paste0(tmp.dir,tasmin.file),write=FALSE)

  pr.scale <- 1
  temp.offset <- 0
  pr.units <- ncatt_get(pr.nc,'pr')$units
  if (pr.units == 'kg m-2 s-1')
    pr.scale <- 86400
  tx.units <- ncatt_get(tasmax.nc,'tasmax')$units
  if (tx.units == 'K')
    temp.scale <- 273
  print('Reading attributes')
  ##Attributes to retain
  lon <- ncvar_get(pr.nc,'lon')
  lat <- ncvar_get(pr.nc,'lat')  
  n.lon <- length(lon)
  n.lat <- length(lat)

  ##Combine the dates
  time.atts <- ncatt_get(pr.nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  
  tasmax.values <- ncvar_get(tasmax.nc,'time')
  tasmax.dates <- as.Date(as.character(past.origin + tasmax.values*86400))

  if (grepl('HadGEM2',tasmax.file)) {
     dx <- which(is.na(tasmax.dates))
     tasmax.dates[dx] <- paste0(format(tasmax.dates[dx-5],'%Y'),'-02-28')
  }

  n.time <- length(tasmax.values)

  print('Latitude loop')

  ##--------------------------------------------------------------
  ##Compute climdex values and load into newly created climdex netcdf
  for (i in 1:n.lat) {
    print(paste('Lat: ',i,' in ',n.lat,sep=''))
    pr.subset <- ncvar_get(pr.nc,'pr',start=c(1,i,1),count=c(n.lon,1,-1))*pr.scale
    pr.list <- separate_into_list(pr.subset)          
    snowdepth.list <- lapply(pr.list,function(x){return(as.numeric(rep(NA,n.time)))})
    swe.list <- lapply(pr.list,function(x){return(as.numeric(rep(NA,n.time)))})
    lat_deg <- lat[i]
    ##Check for NA values
    flag <- is.na(pr.subset[,1000])

    if (sum(flag)==length(flag)) {
      print('All NA values')

    } else {  ##No or some NA values
      if (sum(flag) == 0) {
         print('No NA Values')
      } else {
         print('Some NA Values')
      }

      print('Reading tas data')
      tasmax.subset <- ncvar_get(tasmax.nc,'tasmax',start=c(1,i,1),count=c(n.lon,1,-1))-temp.offset
      tasmin.subset <- ncvar_get(tasmin.nc,'tasmin',start=c(1,i,1),count=c(n.lon,1,-1))-temp.offset      
       
      ##Deal with NA values by omitting cells
      ix <- which(!flag)
      ix.len <- length(ix)
      print(ix.len)

      tasmax.list <- separate_into_list(tasmax.subset[ix,])          
      tasmin.list <- separate_into_list(tasmin.subset[ix,])          
      pr.list <- separate_into_list(pr.subset[ix,])          

      ##-----------------------
      ##Slope and Aspects
      aspect.list <- as.list(rep(0,length(pr.list)))
      slopes.list <- as.list(rep(0,length(pr.list)))
      for (j in 1:ix.len) {
         x <- ix[j]
         aspect.list[[j]] <- bc.aspects[which.min(abs(lon[x]-bc.lon)),which.min(abs(lat[i]-bc.lat))]
         slopes.list[[j]] <- bc.slopes[which.min(abs(lon[x]-bc.lon)),which.min(abs(lat[i]-bc.lat))]
      }
      ##-----------------------
      ##Calibration Parameters
      scale.list <- as.list(rep(0,length(pr.list)))
      slope.list <- as.list(rep(0,length(pr.list)))
      freq.list <- as.list(rep(0,length(pr.list)))

      for (j in 1:ix.len) {
         x <- ix[j]
         scale.list[[j]] <- cal.scales[which.min(abs(lon[x]-cal.lon)),which.min(abs(lat[i]-cal.lat))]
         slope.list[[j]] <- cal.slopes[which.min(abs(lon[x]-cal.lon)),which.min(abs(lat[i]-cal.lat))]
         freq.list[[j]] <- cal.freqs[which.min(abs(lon[x]-cal.lon)),which.min(abs(lat[i]-cal.lat))]
      }

      lat.list <- as.list(rep(lat[i],length(pr.list)))
      dates.list <-rep(list(as.character(tasmax.dates)),length(pr.list))
        
      print('At snow model')
      ptm <- proc.time()          

      snow.objects <- foreach(precip_mm=pr.list,Tmax_C=tasmax.list,Tmin_C=tasmin.list,Date=dates.list,lat_deg=lat.list,
                              slope=slopes.list,aspect=aspect.list,
                              cal_scale=scale.list,cal_slope=slope.list,cal_freq=freq.list,.inorder=TRUE) %dopar% {
                              snow.final <- snow_melt(precip_mm,Tmax_C,Tmin_C,Date,lat_deg,slope,aspect,cal_scale,cal_slope,cal_freq)
                              }          
      print('Done with snow model')
      print('Elapsed time')
      print(proc.time()-ptm)

      snowdepth.values <- lapply(snow.objects,function(x){return(x$snowdepth)})
      swe.values <- lapply(snow.objects,function(x){return(x$swe)})
                  
      snowdepth.list[!flag] <- snowdepth.values
      swe.list[!flag] <- swe.values

      print('Writing')
      snowdepth.matrix <- matrix(unlist(snowdepth.list),nrow=n.lon,ncol=n.time,byrow=TRUE)
      ncvar_put(depth.nc,varid='snowdepth',vals=snowdepth.matrix,
                  start=c(1,i,1),count=c(-1,1,-1))
      swe.matrix <- matrix(unlist(swe.list),nrow=n.lon,ncol=n.time,byrow=TRUE)
      ncvar_put(swe.nc,varid='swe',vals=swe.matrix,
                  start=c(1,i,1),count=c(-1,1,-1))

    }

  }
  print('Total Elapsed time')
  print(proc.time()-atm)
  nc_close(depth.nc)
  nc_close(swe.nc)

  nc_close(pr.nc)
  nc_close(tasmax.nc)
  nc_close(tasmin.nc)
##  closeCluster(cl)
  ##--------------------------------------------------------------
}

##**************************************************************************************

##gcm <- 'HadGEM2-CC'
##run <- 'r1i1p1'
##tmpdir <- '/local_temp/ssobie/snow/'


args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}
reanalysis <- TRUE

if (reanalysis) {
  ##--------------
  ##For reanalysis 
  snow.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/',gcm,'/')  ##For empty snow files
  data.dir <-  paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/',gcm,'/') ##For GCM-PRISM files    
  write.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/calibrated/' ##Write location   

  pr.file <- paste0('pr_gcm_prism_',gcm,'_19790101-20181031.nc')
  tasmax.file <- paste0('tasmax_gcm_prism_',gcm,'_19790101-20181031.nc')
  tasmin.file <- paste0('tasmin_gcm_prism_',gcm,'_19790101-20181031.nc')

  snow.file <- paste0('snowdepth_BCCAQ2-PRISM_',gcm,'_19790101-20181031.nc')
  swe.file <- paste0('swe_BCCAQ2-PRISM_',gcm,'_19790101-20181031.nc')

} else {

  ##-------------
  ##For GCMS
  data.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/',gcm,'/')
  snow.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/' 
  write.dir <- data.dir

  pr.file <- paste0('pr_gcm_prism_',gcm,'_allBC_TPS_1950-2100.nc')
  tasmax.file <- paste0('tasmax_gcm_prism_',gcm,'_allBC_TPS_1950-2100.nc')
  tasmin.file <- paste0('tasmin_gcm_prism_',gcm,'_allBC_TPS_1950-2100.nc')

  snow.file <- paste0('snowdepth_BCCAQ2-PRISM_',gcm,'_19500101-21001231.nc')
  swe.file <- paste0('swe_BCCAQ2-PRISM_',gcm,'_19500101-21001231.nc')

}

tmp.dir <- paste(tmpdir,'/cal_snow_',gcm,'/')

  if (!file.exists(tmp.dir)) {
    dir.create(tmp.dir,recursive=T)
  }
                       
  file.copy(from=paste0(data.dir,pr.file),to=tmp.dir,overwrite=TRUE)
  file.copy(from=paste0(data.dir,tasmax.file),to=tmp.dir,overwrite=TRUE)
  file.copy(from=paste0(data.dir,tasmin.file),to=tmp.dir,overwrite=TRUE)
  file.copy(from=paste0(snow.dir,snow.file),to=tmp.dir,overwrite=TRUE)
  file.copy(from=paste0(snow.dir,swe.file),to=tmp.dir,overwrite=TRUE)

  second <- snow_for_model(var.names=c('snowdepth','swe'),
                           pr.file,tasmax.file,tasmin.file,snow.file,swe.file,
                           data.dir,write.dir,tmp.dir)
  file.copy(from=paste0(tmp.dir,snow.file),to=write.dir,overwrite=TRUE)
  file.copy(from=paste0(tmp.dir,swe.file),to=write.dir,overwrite=TRUE)

  file.remove(paste0(tmp.dir,pr.file))
  file.remove(paste0(tmp.dir,tasmax.file))
  file.remove(paste0(tmp.dir,tasmin.file))
  file.remove(paste0(tmp.dir,snow.file))
  file.remove(paste0(tmp.dir,swe.file))

##**************************************************************************************



