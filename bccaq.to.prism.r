##******************************************************************************
##******************************************************************************

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

library(ncdf4)
library(PCICt)
library(fields)
##******************************************************************************
################################################################################
##Calculate BCCAQ anomalies


##-----------------------------------------------------------------
retrieve_time_info <- function(nc,interval) {
  dates <- netcdf.calendar(nc)
  ymax <-  gsub('-','',format(head(dates,1),'%Y'))
  ymin <-  gsub('-','',format(tail(dates,1),'%Y'))

  bnds <- strsplit(interval,'-')[[1]]
  yst <- head(grep(bnds[1],dates),1)
  yen <- tail(grep(bnds[2],dates),1)
  yct <- yen-yst+1

  rv <- list(dates=dates,ymax=ymax,ymin=ymin,
             yst=yst,yen=yen,yct=yct)
  return(rv)
}

##-----------------------------------------------------------------

baseline_climatologies <- function(var.name,var.clim,monthly.fac,monthly.ts.fac,mon.facs) {
  var.mon <- apply(var.clim,c(1,2),function(x,fac){tapply(x,fac,mean,na.rm=T)},monthly.fac)
  if (var.name=='pr') {
    var.test <- apply(var.clim,c(1,2),function(x,fac){tapply(x,fac,sum,na.rm=T)},monthly.ts.fac)
    var.test[var.test <=0] <- NA
    var.mon <-  apply(var.test,c(2,3),function(x,fac){tapply(x,fac,mean,na.rm=T)},mon.facs)
    rm(var.test)
  }
  rm(var.clim)
  return(var.mon)
}

##-----------------------------------------------------------------

create_anoms <- function(var.name,nc,anc,var.dates,yst,yen) {

  ###Baseline time factors
  clim.monthly.fac <- as.factor(format(var.dates[yst:yen],'%m'))
  clim.monthly.ts.fac <- as.factor(format(var.dates[yst:yen],'%Y-%m'))
  clim.mon.facs <- as.factor(format(as.Date(paste(levels(clim.monthly.ts.fac),'-01',sep='')),'%m'))

  lat <- ncvar_get(nc,'lat')
  n.lat <- length(lat)
  l.width <- 2
  l.seq <- seq(1,368,by=l.width) ##Specific to PNWNAmet

  ##seq(1,192,by=8) ###For ERA5
  lon <- ncvar_get(nc,'lon')
  n.lon <- length(lon)
  time <- ncvar_get(nc,'time')
  n.time <- length(time)

  monthly.fac <- as.factor(format(var.dates,'%m'))

  for (ltx in l.seq) {
    print(paste0('Subset: ',ltx,' in 368')) ###192'))
    var.data <- ncvar_get(nc,var.name,start=c(1,ltx,1),count=c(-1,l.width,-1))
    if (var.name=='pr') {
       var.data[var.data<0] <- 0
    }
    var.clim <- var.data[,,yst:yen]
    var.mon <- baseline_climatologies(var.name,var.clim,
                                      clim.monthly.fac,
                                      clim.monthly.ts.fac,
                                      clim.mon.facs)
    rm(var.clim)
    ##Load full time series and take anomalies from this
    var.anoms <- var.data*0
    for(mn in 1:12) {
      print(mn)
      var.mean <- var.mon[mn,,]

      var.ix <- which(monthly.fac==sprintf('%02d',mn))
      mlen <- length(var.ix)
      for (i in 1:mlen) {
        ix <- var.ix[i]
        if (var.name=='pr')
          var.anoms[,,ix] <- var.data[,,ix]/var.mean
        if (grepl('tas',var.name))
          var.anoms[,,ix] <- var.data[,,ix] - var.mean
      }
    }
    rm(var.data)
    rm(var.mean)
    ncvar_put(anc,varid=var.name,vals=var.anoms,start=c(1,ltx,1),count=c(n.lon,l.width,n.time))
    rm(var.anoms)
  }
  nc_close(nc)
  nc_close(anc)
  gc()
}

##-----------------------------------------------------------------


bccaq_anomalies <- function(var.name,gcm,scenario,interval,base.dir,tmp.dir) {
  print(paste('BCCAQ Anomalies: ',gcm,', ',var.name,sep=''))

  gcm.dir <- paste(base.dir,gcm,'/',sep='')
  var.file <- list.files(path=gcm.dir,pattern=paste(var.name,'_BCCAQ2_',sep=''))
  if (!grepl('TPS',var.file)) {
     stop('Specific to PNWNAmet. Make sure the 368 loop matches the dimensions of the file')
  }
  print(var.file)
  file.copy(from=paste0(base.dir,gcm,'/',var.file),to=tmp.dir,overwrite=TRUE)
  Sys.sleep(5)
  var.tmp.file <- paste0(tmp.dir,'/',var.file)

  ##For mean subset 1981-2010
  anoms.file <- gsub(pattern='_BCCAQ2_',replacement='_anoms_BCCAQ2_',var.file)
  anoms.tmp.file <- paste0(tmp.dir,'/',anoms.file)
  print('Anoms file')
  print(anoms.file)
  file.copy(from=var.tmp.file,to=anoms.tmp.file,overwrite=T)
  Sys.sleep(5)
  if (!file.exists(anoms.tmp.file))
    stop('Anomaly copy failed')

  nc <- nc_open(var.tmp.file)
  anc <- nc_open(anoms.tmp.file,write=TRUE)
  print(anc)
  ##Time series and climatology subset
  time.info <- retrieve_time_info(nc,interval)
  print('Calculate anomalies')
  anoms.past <- create_anoms(var.name,nc,anc,
                             time.info$dates,
                             time.info$yst,
                             time.info$yen)
  return(anoms.file)
}

##-----------------------------------------------------------------

interp_bccaq <- function(var.name,gcm,interval,tmp.dir,anoms.file,grid.file) {
  print(paste('Interpolate Anomalies: ',gcm,', ',var.name,', ',interval,sep=''))
  interp.file <- gsub(pattern='_anoms_BCCAQ2_',replacement='_anoms_interp_',anoms.file)
  system(paste('cdo -s remapbil,',grid.file,' ',tmp.dir,anoms.file,' ',tmp.dir,interp.file,sep=''))
  return(interp.file)
}

##-----------------------------------------------------------------

daily_prism_scale <- function(var.name,gcm,tmp.dir,interp.file,prism.dir) {
  
  prism.var <- switch(var.name,
                      pr='pr',
                      tasmax='tmax',
                      tasmin='tmin')
  prism.file <- paste(prism.dir,prism.var,'_monClim_PRISM_MODIS_GRID_198101-201012.nc',sep='')

  adjusted.file <- gsub(pattern='_anoms_interp_',replacement='_gcm_prism_',interp.file)
  file.copy(from=paste0(tmp.dir,interp.file),to=paste0(tmp.dir,adjusted.file),overwrite=T)  
  ##PRISM climatologies
  pnc <- nc_open(prism.file)
  prism.clim <- ncvar_get(pnc,prism.var)
  nc_close(pnc)
  
  bnc <- nc_open(paste0(tmp.dir,interp.file))
  anc <- nc_open(paste0(tmp.dir,adjusted.file),write=TRUE)
  
  time.atts <- ncatt_get(bnc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(bnc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  var.dates <- origin.pcict + time.values*86400

  monthly.fac <- as.factor(format(var.dates,'%m'))

  ##Break up the data into pieces to be manageable
  tlen <- length(time.values)
  sqnce <- seq(0,100000,by=4000)
  sx <- findInterval(tlen,sqnce)
  itvls <- c(sqnce[1:sx],tlen-sqnce[sx]+sqnce[sx])

  nlon <- bnc$dim$lon$len
  nlat <- bnc$dim$lat$len
  addon <- matrix(NA,nrow=nlon,ncol=9)
  
  for (s in 1:sx) {
    st <- itvls[s] + 1
    en <- itvls[s+1]
    len <- en-st+1
    print(paste(st,' to ',en,' of ',tlen,sep=''))
    var.data <- ncvar_get(bnc,var.name,start=c(1,1,st),count=c(-1,-1,len))

    var.adjust <- var.data*0  
    fac.sub <- monthly.fac[st:en]
    for(mn in 1:12) {
      print(mn)
      prism.mean <- prism.clim[,,mn] ##cbind(addon,prism.clim[,,mn])
      var.ix <- which(sprintf('%02d',fac.sub)==sprintf('%02d',mn))
      mlen <- length(var.ix)
      for (i in 1:mlen) {
        ix <- var.ix[i]
        var.sub <- var.adjust[,,ix]
        if (var.name=='pr') {
          var.sub <- var.data[,,ix]*prism.mean
        }
        if (grepl('tas',var.name)) {
          var.sub <- var.data[,,ix] + prism.mean
        }
        var.adjust[,,ix] <- var.sub
      }##Loop over indices    
    }##Loop over Months
    ##Flag the NA values
    ##browser()
    ##if (var.name=='pr') {
    ##  var.adjust[flags] <- NA
    ##}
    ##var.adjust[,1:9,] <- NA
    ncvar_put(anc,varid=var.name,vals=var.adjust,start=c(1,1,st),count=c(nlon,nlat,len))
  }##Loop over data pieces
  rm(var.adjust)
  nc_close(bnc)
  nc_close(anc)
  return(adjusted.file)
}
################################################################################
##******************************************************************************

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

###gcm <- 'ACCESS1-0'
###var.name <- 'tasmax'
interval <- '1981-2010' ##Baseline for the anomalies
tmp.dir <- tmpdir
var.name <- varname
###tmp.dir <- '/local_temp/ssobie/prism/'

if (!file.exists(tmp.dir))
   dir.create(tmp.dir,recursive=T)

base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/'
prism.dir <- '/storage/data/projects/rci/data/winter_sports/PRISM/'
grid.file <- '/storage/home/ssobie/grid_files/van.whistler.modis.grid.txt'

anoms.file <- bccaq_anomalies(var.name,gcm,scenario,interval,base.dir,tmp.dir)
##print('Copy anomalies back')
##file.copy(from=paste0(tmp.dir,'/',anoms.file),to=paste0(base.dir,gcm,'/'),overwrite=T)
##Sys.sleep(5)

##interp.file <- interp_bccaq(var.name,gcm,'1950-2100',tmp.dir,anoms.file,grid.file)
interp.file <- interp_bccaq(var.name,gcm,'1945-2012',tmp.dir,anoms.file,grid.file)
##file.copy(from=paste0(tmp.dir,interp.file),to=paste0(base.dir,gcm,'/'),overwrite=TRUE)
##Sys.sleep(5)

adjusted.file <- daily_prism_scale(var.name,gcm,tmp.dir,interp.file,prism.dir)
file.copy(from=paste0(tmp.dir,adjusted.file),to=paste0(base.dir,gcm,'/'),overwrite=TRUE)


file.remove(from=paste0(tmp.dir,anoms.file))
file.remove(from=paste0(tmp.dir,interp.file))
file.remove(from=paste0(tmp.dir,adjusted.file))








