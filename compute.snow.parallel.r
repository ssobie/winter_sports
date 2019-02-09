##Script to calculate the climdex indices from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file

##Updated version from compute.climdex.bccaq.r
##This computes all the climdex variables


##source('/storage/data/projects/rci/assessments/code/snow.model.r',chdir=T)
##source('/storage/home/ssobie/code/repos/winter_sports/test.snow.model.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/snow.model.functional.r',chdir=T)

library(ncdf4)
library(PCICt)
library(zoo)

library(doParallel)
registerDoParallel(cores=4) # or some other number that you're comfortable with.
library(foreach)

##---------------------------------------------------------------
separate.into.list <- function(data.subset) {
    return(lapply(seq_len(nrow(data.subset)), function(k) data.subset[k,]))
}


snow.for.model <- function(var.names,
                           pr.file,tasmax.file,tasmin.file,snow.file,swe.file,
                           data.dir,write.dir,tmp.dir) {

  atm <- proc.time()

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

  clim.files <- list(paste0(tmp.dir,snow.file),
                     paste0(tmp.dir,swe.file))
  print(clim.files)
  clim.ncs <- lapply(clim.files,nc_open,write=TRUE)

  ##--------------------------------------------------------------
  print('Reading past')
  pr.past.nc <- nc_open(paste0(tmp.dir,pr.file),write=FALSE)
  tasmax.past.nc <- nc_open(paste0(tmp.dir,tasmax.file),write=FALSE)
  tasmin.past.nc <- nc_open(paste0(tmp.dir,tasmin.file),write=FALSE)

  pr.scale <- 1
  temp.offset <- 0
  pr.units <- ncatt_get(pr.past.nc,'pr')$units
  if (pr.units == 'kg m-2 s-1')
    pr.scale <- 86400
  tx.units <- ncatt_get(tasmax.past.nc,'tasmax')$units
  if (tx.units == 'K')
    temp.scale <- 273
  print('Reading attributes')
  ##Attributes to retain
  lon <- ncvar_get(pr.past.nc,'lon')
  lat <- ncvar_get(pr.past.nc,'lat')  
  n.lon <- length(lon)
  n.lat <- length(lat)

  ##Combine the dates
  time.atts <- ncatt_get(pr.past.nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  
  pr.past.values <- ncvar_get(pr.past.nc,'time')
  pr.time.values <- pr.past.values
  pr.dates <- past.origin + pr.past.values*86400

  tasmax.past.values <- ncvar_get(tasmax.past.nc,'time')
  tasmax.time.values <- tasmax.past.values
  tasmax.dates <- as.Date(as.character(past.origin + tasmax.past.values*86400))

  tasmin.past.values <- ncvar_get(tasmin.past.nc,'time')
  tasmin.time.values <- tasmin.past.values

  print('Latitude loop')

  ##--------------------------------------------------------------
  ##Compute climdex values and load into newly created climdex netcdf
  for (i in 1:n.lat) {
    print(paste('Lat: ',i,' in ',n.lat,sep=''))
    pr.past.subset <- ncvar_get(pr.past.nc,'pr',start=c(1,i,1),count=c(n.lon,1,-1))*pr.scale
    lat_deg <- lat[i]
    ##Check for NA values
    pr.na <- is.na(pr.past.subset[,1])
    if (sum(pr.na)==length(pr.na)) {
      print('All NA values')
      for (k in seq_along(var.names)) {
        ncol <- length(tasmax.dates)
        snow.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)        
      }
    } else {
      print('Reading tas data')
      tasmax.past.subset <- ncvar_get(tasmax.past.nc,'tasmax',start=c(1,i,1),count=c(n.lon,1,-1))-temp.offset
      tasmin.past.subset <- ncvar_get(tasmin.past.nc,'tasmin',start=c(1,i,1),count=c(n.lon,1,-1))-temp.offset
      
      pr.subset <- pr.past.subset
      tasmax.subset <- tasmax.past.subset
      tasmin.subset <- tasmin.past.subset

      snowdepth.list  <- vector(mode='list',length=n.lon)
      swe.list  <- vector(mode='list',length=n.lon)
      
      ##Deal with NA values by omitting cells
      if (sum(pr.na)>0) {
        print('Some NA values present')
          ix <- which(!pr.na)
          ix.len <- length(ix)
        print(ix.len)

          tasmax.list <- separate.into.list(tasmax.subset[ix,])          
          tasmin.list <- separate.into.list(tasmin.subset[ix,])          
          pr.list <- separate.into.list(pr.subset[ix,])          

          aspect.list <- as.list(rep(0,length(pr.list)))
          slopes.list <- as.list(rep(0,length(pr.list)))
          for (j in 1:ix.len) {
            x <- ix[j]
            aspect.list[[j]] <- bc.aspects[which.min(abs(lon[x]-bc.lon)),which.min(abs(lat[i]-bc.lat))]
            slopes.list[[j]] <- bc.slopes[which.min(abs(lon[x]-bc.lon)),which.min(abs(lat[i]-bc.lat))]
          }
          lat.list <- as.list(rep(lat[i],length(pr.list)))
          dates.list <-rep(list(as.character(tasmax.dates)),length(pr.list))
        
          print('At snow model')
          ptm <- proc.time()          
if (1==1) {
          snow.objects <- foreach(i=1:n.lon,
                                  precip_mm=pr.list,
                                  Tmax_C=tasmax.list,
                                  Tmin_C=tasmin.list,
                                  Date=dates.list,
                                  lat_deg=lat.list,
                                  slope=slopes.list,
                                  aspect=aspect.list,
                                  .inorder=TRUE) %dopar% {
                                  ##.options.mpi=mpi.options) %dopar% {
                                    ##.export=c('tasmax.dates','lat_deg')
                                    snow.final <- snow.melt(precip_mm,Tmax_C,Tmin_C,Date,lat_deg,slope,aspect)                          
                                  }
}           

##          snow.objects <- mapply(snow.melt,pr.list,tasmax.list,tasmin.list,
##                                 MoreArgs=list(Date=tasmax.dates,lat=lat[i]),SIMPLIFY=FALSE)
          print('Done with snow model')
          print('Elapsed time')
          print(proc.time()-ptm)
                  
          nx <- which(pr.na)
          for (n in 1:length(nx)) {
            y <- nx[n]
            snowdepth.list[[y]] <- rep(NA,length(tasmax.dates))
            swe.list[[y]] <- rep(NA,length(tasmax.dates))
          }

          for (j in 1:ix.len) {
            x <- ix[j]
            snow.data <- snow.objects[[j]]
            snowdepth.list[[x]] <- snow.data$snowdepth
            swe.list[[x]] <- snow.data$swe
          }
          snow.list <- list(snowdepth=snowdepth.list,
                            swe=swe.list)

        } else {
          print('No NA values present')
          
          tasmax.list <- separate.into.list(tasmax.subset[ix,])          
          tasmin.list <- separate.into.list(tasmin.subset[ix,])          
          pr.list <- separate.into.list(pr.subset[ix,])          

          aspect.list <- as.list(rep(lat[i],length(pr.list)))
          slopes.list <- as.list(rep(lat[i],length(pr.list)))
          lat.list <- as.list(rep(lat[i],length(pr.list)))
          dates.list <-rep(list(as.character(tasmax.dates)),length(pr.list))


          for (j in 1:n.lon) {
            aspect.list[[j]] <- bc.aspects[which.min(abs(lon[j]-bc.lon)),which.min(abs(lat[i]-bc.lat))]
            slopes.list[[j]] <- bc.slopes[which.min(abs(lon[j]-bc.lon)),which.min(abs(lat[i]-bc.lat))]
          }

          print('At snow model')
          ptm <- proc.time()
if (1==1) {
          snow.objects <- foreach(i=1:n.lon,          
                                  precip_mm=pr.list,
                                  Tmax_C=tasmax.list,
                                  Tmin_C=tasmin.list,
                                  Date=dates.list,
                                  lat_deg=lat.list,
                                  slope=slopes.list,
                                  aspect=aspect.list,
                                  .inorder=TRUE) %dopar% {
                                  ##.options.mpi=mpi.options) %dopar% {
                                    snow.final <- snow.melt(precip_mm,Tmax_C,Tmin_C,Date,lat_deg,slope,aspect)                          
                                  }
}

##          snow.objects <- mapply(snow.melt,pr.list,tasmax.list,tasmin.list,
##                                 MoreArgs=list(Date=tasmax.dates,lat=lat[i]),SIMPLIFY=FALSE)

          print('Done with snow model')
          print('Elapsed time')
          print(proc.time()-ptm)
          ltm <- proc.time()
          for (n in 1:n.lon) {
            snow.data <- snow.objects[[n]]
            snowdepth.list[[n]] <- snow.data$snowdepth
            swe.list[[n]] <- snow.data$swe
          }
          snow.list <- list(snowdepth=snowdepth.list,
                            swe=swe.list)
          print('Loop timing')
          print(proc.time()-ltm)
        }
      print('Writing')
      for (k in seq_along(var.names)) {
        snow.values <- snow.list[[var.names[k]]]
        ncol <- length(tasmax.dates)
        snow.matrix <- matrix(unlist(snow.values),nrow=n.lon,ncol=ncol,byrow=TRUE)

        ncvar_put(clim.ncs[[k]],varid=var.names[k],vals=snow.matrix,
                  start=c(1,i,1),count=c(-1,1,-1))
      }
    }

  }
  print('Total Elapsed time')
  print(proc.time()-atm)
  lapply(clim.ncs,nc_close)

  nc_close(pr.past.nc)
  nc_close(tasmax.past.nc)
  nc_close(tasmin.past.nc)
##  closeCluster(cl)
  ##--------------------------------------------------------------
}

##**************************************************************************************


args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

 
  data.dir <-  paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/',gcm,'/') 
  write.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/' 
  tmp.dir <- tmpdir ##'/local_temp/ssobie/snow/'

  if (!file.exists(tmp.dir)) {
    dir.create(tmp.dir,recursive=T)
  }

  ##gcm <- 'ERA'

  pr.file <- paste0('pr_gcm_prism_',gcm,'_19790101-20181031.nc')
  tasmax.file <- paste0('tasmax_gcm_prism_',gcm,'_19790101-20181031.nc')
  tasmin.file <- paste0('tasmin_gcm_prism_',gcm,'_19790101-20181031.nc')
  snow.file <- paste0('snowdepth_BCCAQ2-PRISM_',gcm,'_19790101-20181031.nc')
  swe.file <- paste0('swe_BCCAQ2-PRISM_',gcm,'_19790101-20181031.nc')

  
  file.copy(from=paste0(data.dir,pr.file),to=tmp.dir,overwrite=TRUE)
  file.copy(from=paste0(data.dir,tasmax.file),to=tmp.dir,overwrite=TRUE)
  file.copy(from=paste0(data.dir,tasmin.file),to=tmp.dir,overwrite=TRUE)
  file.copy(from=paste0(data.dir,snow.file),to=tmp.dir,overwrite=TRUE)
  file.copy(from=paste0(data.dir,swe.file),to=tmp.dir,overwrite=TRUE)

  second <- snow.for.model(var.names=c('snowdepth','swe'),
                           pr.file,tasmax.file,tasmin.file,snow.file,swe.file,
                           data.dir,write.dir,tmp.dir)
  file.copy(from=paste0(tmp.dir,snow.file),to=paste0(write.dir,'snow/'),overwrite=TRUE)
  file.copy(from=paste0(tmp.dir,swe.file),to=paste0(write.dir,'snow/'),overwrite=TRUE)

  file.remove(paste0(tmp.dir,pr.file))
  file.remove(paste0(tmp.dir,tasmax.file))
  file.remove(paste0(tmp.dir,tasmin.file))
  file.remove(paste0(tmp.dir,snow.file))
  file.remove(paste0(tmp.dir,swe.file))

##**************************************************************************************



