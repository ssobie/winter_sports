##******************************************************************************
##******************************************************************************

library(ncdf4)
library(PCICt)
library(fields)
##******************************************************************************
################################################################################
##Calculate BCCAQ anomalies

create.anoms <- function(var.name,file,var.mon) {

  ##For mean subset 1971-2000
  anoms.file <- gsub(pattern='_day_',replacement='_anoms_',file)
  file.copy(from=file,to=anoms.file,overwrite=T)
  Sys.sleep(5)
  if (!file.exists(anoms.file))
    browser()
  
  nc <- nc_open(file)
  var.data <- ncvar_get(nc,var.name)
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  var.dates <- origin.pcict + time.values*86400/24
  monthly.fac <- as.factor(format(var.dates,'%m'))
  
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
  nc_close(nc)

  ##Fix the edges for interpolations
  ncol <- dim(var.data)[2]
  for (j in 1:(ncol-1)) {
    ix <- is.na(var.anoms[,j,1])
    var.anoms[ix,j,] <- var.anoms[ix,(j+1),]    
  }
  
  if (var.name=='pr')
    var.anoms[is.na(var.anoms)] <- 0


  anc <- nc_open(anoms.file,write=TRUE)
  ncvar_put(anc,varid=var.name,vals=var.anoms)
  nc_close(anc)
  gc()
}

bccaq.anomalies <- function(var.name,gcm,base.dir) {
  print(paste('BCCAQ Anomalies: ',gcm,', ',var.name,sep=''))

  base.file <- list.files(path=paste(base.dir,'baseline/',gcm,sep=''),pattern=paste(var.name,'_day_',sep=''),full.name=TRUE)

  gcm.dir <- paste(base.dir,gcm,sep='')
  print(base.file)
  gnc <- nc_open(base.file)
  time.atts <- ncatt_get(gnc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(gnc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  var.dates <- origin.pcict + time.values*86400/24
  monthly.fac <- as.factor(format(var.dates,'%m'))
  monthly.ts.fac <- as.factor(format(var.dates,'%Y-%m'))
  mon.facs <- as.factor(format(as.Date(paste(levels(monthly.ts.fac),'-01',sep='')),'%m'))
  var.data <- ncvar_get(gnc,var.name)
  var.mon <- apply(var.data,c(1,2),function(x,fac){tapply(x,fac,mean,na.rm=T)},monthly.fac)
  if (var.name=='pr') {
    var.data[var.data <=0] <- NA    
    var.test <- apply(var.data,c(1,2),function(x,fac){tapply(x,fac,sum,na.rm=T)},monthly.ts.fac)
    var.mon <-  apply(var.test,c(2,3),function(x,fac){tapply(x,fac,mean,na.rm=T)},mon.facs)
  }

  nc_close(gnc)
  ##------------------------------------  
  var.files <- list.files(path=gcm.dir,pattern=paste(var.name,'_day_',sep=''),full.name=TRUE)
  past.file <- var.files[grep('1979-2016',var.files)]
  ##proj.file <- var.files[grep('2001-2100',var.files)]

  create.anoms(var.name,past.file,var.mon)
  ##create.anoms(var.name,proj.file,var.mon)
}

interp.bccaq <- function(var.name,gcm,interval,base.dir,grid.file) {
  print(paste('Interpolate Anomalies: ',gcm,', ',var.name,', ',interval,sep=''))
  gcm.dir <- paste(base.dir,gcm,sep='')  
  anoms.files <- list.files(path=gcm.dir,pattern=paste(var.name,'_anoms_QDM_',sep=''),full.name=TRUE)
  anoms.file <- anoms.files[grep(interval,anoms.files)]
  interp.file <- gsub(pattern='_anoms_QDM_',replacement='_anoms_interp_',anoms.file)

  system(paste('cdo -s remapbil,',grid.file,' ',anoms.file,' ',interp.file,sep=''))
  gc()  
}

daily.prism.scale <- function(var.name,gcm,interval,base.dir,prism.dir) {
  print(paste('Daily PRISM: ',gcm,', ',var.name,', ',interval,sep=''))
  gcm.dir <- paste(base.dir,gcm,'/',sep='')
  bccaq.files <- list.files(path=gcm.dir,pattern=paste(var.name,'_anoms_interp_',sep=''),full.name=TRUE)
  bccaq.file <- bccaq.files[grep(interval,bccaq.files)]
  prism.var <- switch(var.name,
                      pr='pr',
                      tasmax='tmax',
                      tasmin='tmin')
  prism.file <- paste(prism.dir,prism.var,'_monClim_PRISM_MODIS_GRID_198101-201012.nc',sep='')

  adjusted.file <- gsub(pattern='_anoms_interp_',replacement='_gcm_prism_',bccaq.file)
  file.copy(from=bccaq.file,to=adjusted.file,overwrite=T)  
  ##PRISM climatologies
  pnc <- nc_open(prism.file)
  prism.clim <- ncvar_get(pnc,prism.var)
  nc_close(pnc)
  
  bnc <- nc_open(bccaq.file)
  anc <- nc_open(adjusted.file,write=TRUE)
  
  time.atts <- ncatt_get(bnc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(bnc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  var.dates <- origin.pcict + time.values*86400/24
  monthly.fac <- as.factor(format(var.dates,'%m'))

##Flags for the SW corner
  flag.file <- list.files(path=paste(base.dir,gcm,sep=''),pattern='tasmax_gcm_prism_',full.name=TRUE)[1]
  fnc <- nc_open(flag.file,write=TRUE)
  flag.data <- ncvar_get(fnc,'tasmax',start=c(1,1,1),count=c(-1,-1,1))
  flags <- is.na(flag.data)
  nc_close(fnc)
  
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
      var.ix <- which(fac.sub==sprintf('%02d',mn))
      mlen <- length(var.ix)
      for (i in 1:mlen) {
        ix <- var.ix[i]
        var.sub <- var.adjust[,,ix]
        if (var.name=='pr') {
          var.sub <- var.data[,,ix]*prism.mean
        }
          ##var.adjust[,10:nlat,ix] <- var.data[,10:nlat,ix]*prism.mean
        if (grepl('tas',var.name)) {
          var.sub <- var.data[,,ix] + prism.mean
        }
          ##var.adjust[,10:nlat,ix] <- var.data[,10:nlat,ix] + prism.mean
        var.sub[flags] <- NA        
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
  gc()
}
################################################################################
##******************************************************************************

run.adjust <- function() {

  ##Requires 1971-2000 (or equivalent base period) in the /baseline directory to create anomalies from the full period
  ## (usually 1950-2000 and 2001-2100). Both of these are created using extract.bccaq.gcm.r
  ##Also need the PRISM climatologies (also using extract.bccaq.gcm.r).

  base.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/'
  prism.dir <- '/storage/data/projects/rci/data/winter_sports/PRISM/'
  grid.file <- '/storage/home/ssobie/grid_files/van.whistler.modis.grid.txt'

  var.list <- 'pr' ##c('tasmax','tasmin','pr')

  gcm.list <- c('ACCESS1-0',
                'CCSM4',
                'CanESM2',
                'CNRM-CM5',
                'CSIRO-Mk3-6-0',
                'GFDL-ESM2G',
                'HadGEM2-CC',
                'HadGEM2-ES',  
                'inmcm4',
                'MIROC5',
                'MPI-ESM-LR',
                'MRI-CGCM3')
  gcm.list <- 'NCEP2'

  for (var.name in var.list) {
    print(var.name)
    for (gcm in gcm.list) {
      print(gcm)
      bccaq.anomalies(var.name,gcm,base.dir)
      interp.bccaq(var.name,gcm,'1979-2016',base.dir,grid.file)
         ##interp.bccaq(var.name,gcm,'2001-2100',base.dir,grid.file)
      daily.prism.scale(var.name,gcm,'1979-2016',base.dir,prism.dir)
         ##daily.prism.scale(var.name,gcm,'2001-2100',base.dir)
    }
  }  
}

run.adjust()



