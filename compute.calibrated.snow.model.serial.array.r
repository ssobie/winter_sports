##Script to calculate the climdex indices from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file

##Rprof('snow_model_profile.out')

##This version uses split apart GCM-PRISM and snow files

source('/storage/home/ssobie/code/repos/winter_sports/snow.model.calibrated.r',chdir=T)

library(ncdf4)
library(PCICt)
library(zoo)

##---------------------------------------------------------------
separate_into_list <- function(data.subset) {
    return(lapply(seq_len(nrow(data.subset)), function(k) data.subset[k,]))
}


snow_for_model <- function(var.names,model,
                           pr.file,tasmax.file,tasmin.file,snow.file,swe.file,
                           scale.file,slope.file,freq.file,
                           tmp.dir) {

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
  scale.nc <- nc_open(paste0(tmp.dir,scale.file))
  cal.scales <- ncvar_get(scale.nc,'scale') * -1 ##to correct for kriging with positive values
  cal.lon <- ncvar_get(scale.nc,'lon')
  cal.lat <- ncvar_get(scale.nc,'lat')
  nc_close(scale.nc)

  slope.nc <- nc_open(paste0(tmp.dir,slope.file))
  cal.slopes <- ncvar_get(slope.nc,'slope')
  nc_close(slope.nc)

  freq.nc <- nc_open(paste0(tmp.dir,freq.file))
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
  
  pr.values <- ncvar_get(pr.nc,'time')
  pr.dates <- as.Date(as.character(past.origin + pr.values*86400))

  if (grepl('HadGEM2',pr.file)) {
     dx <- which(is.na(pr.dates))
     pr.dates[dx] <- paste0(format(pr.dates[dx-5],'%Y'),'-02-28')
  }

  n.time <- length(pr.values)

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
      if (grepl('HadGEM2',model)) {
         t.end <- 54000
      } else {
         t.end <- -1
      }
      tasmax.subset <- ncvar_get(tasmax.nc,'tasmax',start=c(1,i,1),count=c(n.lon,1,t.end))-temp.offset
      tasmin.subset <- ncvar_get(tasmin.nc,'tasmin',start=c(1,i,1),count=c(n.lon,1,t.end))-temp.offset      
       
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
      dates.list <-rep(list(as.character(pr.dates)),length(pr.list))

      print('At snow model')
      ptm <- proc.time()          
      snow.objects <- vector(mode='list',length=ix.len)
      for (j in 1:ix.len) {
         snow.objects[[j]] <- snow_melt(precip_mm=pr.list[[j]],Tmax_C=tasmax.list[[j]],Tmin_C=tasmin.list[[j]],
                                        Date=dates.list[[j]],lat_deg=lat.list[[j]],
                                        cal_scale=scale.list[[j]],cal_slope=slope.list[[j]],cal_freq=freq.list[[j]],
                                        slope=slopes.list[[j]],aspect=aspect.list[[j]])
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

  ##--------------------------------------------------------------
}

##**************************************************************************************
testing <- FALSE

if (testing) {

   array_value <- 1

   model <- 'HadGEM2-CC'
   reanalysis <- 'PNWNAmet'
   tmpdir <- '/local_temp/ssobie/snowcal_testing/'
   caldir <- "/storage/data/projects/rci/data/winter_sports/"
   datadir <- "/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/"
   snowdir <- "/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/templates/"
   writedir <- "/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/"

   prfile <- list.files(path=paste0(datadir,model,'/pr_',model,'_split10/'),pattern='pr_L')[array_value]
   txfile <- list.files(path=paste0(datadir,model,'/tasmax_',model,'_split10/'),pattern='tasmax_L')[array_value]
   tnfile <- list.files(path=paste0(datadir,model,'/tasmin_',model,'_split10/'),pattern='tasmin_L')[array_value]

   snowfile <- list.files(path=paste0(snowdir,model,'/snowdepth_',model,'_split10/'),pattern='snowdepth_L')[array_value]
   swefile <- list.files(path=paste0(snowdir,model,'/swe_',model,'_split10/'),pattern='swe_L')[array_value]

   scalefile <- list.files(path=paste0(caldir,'scale_',reanalysis,'_split10/'),pattern='scale_L')[array_value]
   slopefile <- list.files(path=paste0(caldir,'slope_',reanalysis,'_split10/'),pattern='slope_L')[array_value]
   freqfile <- list.files(path=paste0(caldir,'freq_',reanalysis,'_split10/'),pattern='freq_L')[array_value]

} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
   }
}


tmp.dir <- paste0(tmpdir,'/cal_snow_',model,'/')
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=T)
}

file.copy(from=paste0(datadir,model,'/pr_',model,'_split10/',prfile),to=tmp.dir,overwrite=TRUE)
file.copy(from=paste0(datadir,model,'/tasmax_',model,'_split10/',txfile),to=tmp.dir,overwrite=TRUE)
file.copy(from=paste0(datadir,model,'/tasmin_',model,'_split10/',tnfile),to=tmp.dir,overwrite=TRUE)

file.copy(from=paste0(snowdir,model,'/snowdepth_',model,'_split10/',snowfile),to=tmp.dir,overwrite=TRUE)
file.copy(from=paste0(snowdir,model,'/swe_',model,'_split10/',swefile),to=tmp.dir,overwrite=TRUE)

file.copy(from=paste0(caldir,'scale_',reanalysis,'_split10/',scalefile),to=tmp.dir,overwrite=TRUE)
file.copy(from=paste0(caldir,'slope_',reanalysis,'_split10/',slopefile),to=tmp.dir,overwrite=TRUE)
file.copy(from=paste0(caldir,'freq_',reanalysis,'_split10/',freqfile),to=tmp.dir,overwrite=TRUE)

 
second <- snow_for_model(var.names=c('snowdepth','swe'),model,
                         prfile,txfile,tnfile,snowfile,swefile,
                         scalefile,slopefile,freqfile,
                         tmp.dir)
copy.dir <- paste0(writedir,'calibrated_',model,'_',reanalysis,'_prism_tps/')
if (!file.exists(copy.dir)) {
   dir.create(copy.dir,recursive=T)
}

swe.copy.dir <- paste0(copy.dir,'swe_',model,'_split10/')
if (!file.exists(swe.copy.dir)) {
   dir.create(swe.copy.dir,recursive=T)
}

snow.copy.dir <- paste0(copy.dir,'snowdepth_',model,'_split10/')
if (!file.exists(snow.copy.dir)) {
   dir.create(snow.copy.dir,recursive=T)
}


file.copy(from=paste0(tmp.dir,snowfile),to=snow.copy.dir,overwrite=TRUE)
file.copy(from=paste0(tmp.dir,swefile),to=swe.copy.dir,overwrite=TRUE)

file.remove(paste0(tmp.dir,prfile))
file.remove(paste0(tmp.dir,txfile))
file.remove(paste0(tmp.dir,tnfile))
file.remove(paste0(tmp.dir,snowfile))
file.remove(paste0(tmp.dir,swefile))
file.remove(paste0(tmp.dir,scalefile))
file.remove(paste0(tmp.dir,slopefile))
file.remove(paste0(tmp.dir,freqfile))

##Rprof(NULL)

                       

##**************************************************************************************



