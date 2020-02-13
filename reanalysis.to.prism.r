##******************************************************************************
##******************************************************************************

library(ncdf4)
library(PCICt)
library(fields)
##******************************************************************************
################################################################################
##Calculate BCCAQ anomalies

create.anoms <- function(var.name,tmp.dir,file,var.mon) {

  ##For mean subset 1971-2000
  anoms.file <- gsub(pattern='_day_',replacement='_anoms_',file)
  file.copy(from=paste0(tmp.dir,file),to=paste0(tmp.dir,anoms.file),overwrite=T)
  Sys.sleep(5)
  if (!file.exists(paste0(tmp.dir,anoms.file)))
    browser()
  
  nc <- nc_open(paste0(tmp.dir,file))
  var.data <- ncvar_get(nc,var.name)
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  

  if (grepl('NCEP2',file)) {
     var.dates <- origin.pcict + time.values*3600                            
  }
  if (grepl('ERA',file)) {
     var.dates <- origin.pcict + time.values*86400
  }

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


  anc <- nc_open(paste0(tmp.dir,anoms.file),write=TRUE)
  ncvar_put(anc,varid=var.name,vals=var.anoms)
  nc_close(anc)
  return(anoms.file)
}

bccaq.anomalies <- function(var.name,gcm,tmp.dir,base.file,past.file) {
  print(paste('BCCAQ Anomalies: ',gcm,', ',var.name,sep=''))

  print(base.file)
  gnc <- nc_open(paste0(tmp.dir,base.file))
  time.atts <- ncatt_get(gnc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(gnc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  if (grepl('NCEP2',base.file)) {
     var.dates <- origin.pcict + time.values*3600                            
  }
  if (grepl('ERA',base.file)) {
     var.dates <- origin.pcict + time.values*86400
  }

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
##browser()
  nc_close(gnc)
  ##------------------------------------  

  anoms.file <- create.anoms(var.name,tmp.dir,past.file,var.mon)
  return(anoms.file)
}

interp.bccaq <- function(var.name,gcm,interval,tmp.dir,anoms.file,grid.file) {
  print(paste('Interpolate Anomalies: ',gcm,', ',var.name,', ',interval,sep=''))
  interp.file <- gsub(pattern='_anoms_QDM_',replacement='_anoms_interp_',anoms.file)
  system(paste('cdo -s remapbil,',grid.file,' ',tmp.dir,anoms.file,' ',tmp.dir,interp.file,sep=''))
  return(interp.file)
}

daily.prism.scale <- function(var.name,gcm,tmp.dir,interp.file,prism.dir) {
  print(paste('Daily PRISM: ',gcm,', ',var.name,'19790101-20181031',sep=''))
  
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

  if (grepl('NCEP2',interp.file)) {
     var.dates <- origin.pcict + time.values*3600                            
  }
  if (grepl('ERA',interp.file)) {
     var.dates <- origin.pcict + time.values*86400
  }

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
          ##var.adjust[,10:nlat,ix] <- var.data[,10:nlat,ix]*prism.mean
        if (grepl('tas',var.name)) {
          var.sub <- var.data[,,ix] + prism.mean
        }
          ##var.adjust[,10:nlat,ix] <- var.data[,10:nlat,ix] + prism.mean
        ##var.sub[flags] <- NA        
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

run.adjust <- function() {

  ##Requires 1971-2000 (or equivalent base period) in the /baseline directory to create anomalies from the full period
  ## (usually 1950-2000 and 2001-2100). Both of these are created using extract.bccaq.gcm.r
  ##Also need the PRISM climatologies (also using extract.bccaq.gcm.r).

  tmp.dir <- '/local_temp/ssobie/prism/'
  if (!file.exists(tmp.dir))
     dir.create(tmp.dir,recursive=T)

 base.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/'
  prism.dir <- '/storage/data/projects/rci/data/winter_sports/PRISM/'
  grid.file <- '/storage/home/ssobie/grid_files/van.whistler.modis.grid.txt'

  var.list <- c('tasmax','tasmin','pr')
  gcm.list <- 'ERA'
  
  for (var.name in var.list) {
    print(var.name)
    for (gcm in gcm.list) {
      print(gcm)

      base.file <- paste0(var.name,'_day_QDM_',gcm,'_1981-2010.nc')
      file.copy(from=paste0(base.dir,'baseline/',base.file),to=tmp.dir,overwrite=TRUE)
      past.file <- paste0(var.name,'_day_QDM_',gcm,'_19790101-20181031.nc')
      file.copy(from=paste0(base.dir,gcm,'/',past.file),to=tmp.dir,overwrite=TRUE)
      anoms.file <- bccaq.anomalies(var.name,gcm,tmp.dir,base.file,past.file)
      file.copy(from=paste0(tmp.dir,anoms.file),to=paste0(base.dir,gcm,'/'),overwrite=TRUE)      
      interp.file <- interp.bccaq(var.name,gcm,'1979-2018',tmp.dir,anoms.file,grid.file)
      file.copy(from=paste0(tmp.dir,interp.file),to=paste0(base.dir,gcm,'/'),overwrite=TRUE)
      adjusted.file <- daily.prism.scale(var.name,gcm,tmp.dir,interp.file,prism.dir)
      file.copy(from=paste0(tmp.dir,adjusted.file),to=paste0(base.dir,gcm,'/'),overwrite=TRUE)

    }
  }  
}

run.adjust()



