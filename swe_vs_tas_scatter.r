##Script to plot time series of frost free days for Vancouver Intl.

library(ncdf4)
library(PCICt)
library(rgdal)
library(rgeos)
library(zoo)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

make_swe_vs_temp_anomaly_plot <- function(swe,tas.anoms,
                                          tas.axis,swe.axis,xlim,ylim,
                                          elev,x.axis=FALSE,y.axis=FALSE) {

   plot(tas.anoms[,1],swe.series[[1]]*1000,xlim=xlim,ylim=ylim,pch=16,col='darkgray',
        xlab='',ylab='',axes=FALSE,yaxs='i',cex=1.25)
   abline(h=swe.axis,lwd=0.8,lty=2,col='gray94')
   for (j in seq_along(swe)) {
      points(tas.anoms[,j],swe.series[[j]]*1000,pch=16,col='darkgray',cex=1.25)
   }
   if (x.axis) {
     axis(1,at=tas.axis,label=tas.axis,cex.axis=2.5,mgp=c(4,2,0))
   }
   if (y.axis) {
     axis(2,at=swe.axis,label=swe.axis,cex.axis=2.5,mgp=c(4,2,0))
   }

   boxplot(x=swe.means,at=0,col='blue',add=T,axes=F,cex=2)
   boxplot(x=swe.one,at=1,col='goldenrod',add=T,axes=F,cex=2)
   boxplot(x=swe.two,at=2,col='orange',add=T,axes=F,cex=2)
   boxplot(x=swe.three,at=3,col='red',add=T,axes=F,cex=2)
   text(x=7,y=0.9*ylim[2],paste0('Elevation\n',elev,'m'),cex=2.5)
   box(which='plot')

}

##---------------------------------------------------------------------

prism_elevation <- function(lon.ix,lat.ix,prism.nc) {
   dem <- ncvar_get(prism.nc,'elevation',start=c(lon.ix,lat.ix,1),count=c(1,1,1))
   return(dem)
}

##---------------------------------------------------------------------

read_cell_clim <- function(gcm,var.name,var.file,read.dir,lon.ix,lat.ix,seas.ix) { ##,lonc,latc) {

  ##Coordinate indices
  nc <- nc_open(paste0(read.dir,var.file))
  var.clim <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,seas.ix),count=c(1,1,1))
  var.time <- format(netcdf.calendar(nc),'%Y-%m-%d')

  nc_close(nc)
  return(var.clim)
}

##---------------------------------------------------------------------

read_cell_series <- function(gcm,var.name,var.file,read.dir,lon.ix,lat.ix,seas) { ##,lonc,latc) {

  ##Coordinate indices
  nc <- nc_open(paste0(read.dir,var.file))
  var.series <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  var.time <- format(netcdf.calendar(nc),'%Y-%m-%d')

  seas.ix <- switch(seas,
                    winter='-01-',
                    spring='-04-',
                    summer='-07-',
                    fall='-10-',
                    annual='*')
  time.ix <- grep(seas.ix,var.time,'%m')
  var.sub <- var.series[time.ix]
  time.sub <- var.time[time.ix]

  if (grepl('HadGEM',var.file)) {
     years <- 1950:2100
     had.years <- as.numeric(format(as.Date(time.sub),'%Y'))
     had.ix <- years %in% had.years
     had.fill <- rep(NA,length(years))
     had.fill[had.ix] <- var.sub
     had.time <- paste0(years,'-',format(as.Date(time.sub[20]),'%m-%d'))
     var.sub <- had.fill
     time.sub <- had.time
  }

  nc_close(nc)
  rv <- list(series=var.sub,time=time.sub)
  return(rv)
}

##---------------------------------------------------------------------

base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/'

rcp85.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
                'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')

##Whistler Village
lonc <- -122.954
latc <- 50.1136

##Whistler Mountain
lonc <- -122.950
latc <- 50.0593

##Whistler slope indices
lon.ix <- 125
lat.ix <- 210:217

whistler.elevs <- lat.ix*0

prism.file <- '/storage/data/projects/rci/data/prism/van_whistler_prism_dem_elevations.nc'
prism.nc <- nc_open(prism.file)

for (i in seq_along(lat.ix)) {
   whistler.elevs[i] <- prism_elevation(lon.ix,lat.ix[i],prism.nc)
}
nc_close(prism.nc)

##------------------------------------------------------------
type <- 'swe_seasonal_mean'
seas.ix <- 2 ##Winter
seas <- 'spring'

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/whistler.ski.slope.spring.mean.swe.vs.tas.png')

###png(file=plot.file,width=6,height=8,units='in',res=600,pointsize=6,bg='white')

xlim=c(-3,8)
ylim=c(0,2700)

par(mfrow=c(4,2))
par(mar=c(0,0,0,0),oma=c(8,8,5,5))

for (j in seq_along(lat.ix)) { 
   j <- 7
   swe.data <- vector(mode='list',length=length(rcp85.list))
   swe.one <- swe.two <- swe.three <- rep(0,length(rcp85.list))

   tasmax.data <- vector(mode='list',length=length(rcp85.list))
   tasmin.data <- vector(mode='list',length=length(rcp85.list))

   for (i in seq_along(rcp85.list)) {
      gcm <- rcp85.list[i]
      print(gcm)
      swe.dir <- paste0(base.dir,'snow_model/calibrated_',gcm,'_PNWNAmet_prism_tps/')
      swe.file <- list.files(path=swe.dir,pattern=type)
      swe.data[[i]] <- read_cell_series(gcm,var.name='swe',swe.file,swe.dir,lon.ix,lat.ix[j],seas)

      swe.one.dir <- paste0(swe.dir,'temperature_climatologies/')
      swe.one.file <- list.files(path=swe.one.dir,pattern=type)[1]
      swe.one[i] <- read_cell_clim(gcm,var.name='swe',swe.one.file,swe.one.dir,lon.ix,lat.ix[j],seas.ix)*1000

      swe.two.dir <- paste0(swe.dir,'temperature_climatologies/')
      swe.two.file <- list.files(path=swe.two.dir,pattern=type)[2]
      swe.two[i] <- read_cell_clim(gcm,var.name='swe',swe.two.file,swe.two.dir,lon.ix,lat.ix[j],seas.ix)*1000

      swe.three.dir <- paste0(swe.dir,'temperature_climatologies/')
      swe.three.file <- list.files(path=swe.three.dir,pattern=type)[3]
      swe.three[i] <- read_cell_clim(gcm,var.name='swe',swe.three.file,swe.three.dir,lon.ix,lat.ix[j],seas.ix)*1000

      tas.dir <- paste0(base.dir,gcm,'/tas_climatology/')
      tasmax.file <- list.files(path=tas.dir,pattern=paste0('tasmax_seasonal_mean_gcm_prism_',gcm))
      tasmax.data[[i]] <- read_cell_series(gcm,var.name='tasmax',tasmax.file,tas.dir,lon.ix,lat.ix[j],seas)
      tasmin.file <- list.files(path=tas.dir,pattern=paste0('tasmin_seasonal_mean_gcm_prism_',gcm))
      tasmin.data[[i]] <- read_cell_series(gcm,var.name='tasmin',tasmin.file,tas.dir,lon.ix,lat.ix[j],seas)
   }

   swe.series <- lapply(swe.data,function(x){return(x$series)})

   tas.data <- mapply(FUN=function(x,y){(x$series+y$series)/2},tasmax.data,tasmin.data)
   tas.anoms <- apply(tas.data,2,function(x){x-mean(x[31:61])})
   swe.anoms <- lapply(swe.series,function(x){x-mean(x[31:60])})
   swe.means <- unlist(lapply(swe.series,function(x){mean(x[31:60])}))*1000

   ### make_swe_vs_temp_anomaly_plot(swe.series,tas.anoms,xlim,ylim,whistler.elevs[j])
   x.axis <- y.axis <- FALSE
   if (j %in% c(1,3,5,7)) {y.axis<-TRUE}
   if (j %in% c(7,8)) {x.axis<-TRUE}

   if (j %in% c(1,2)) { 
      ylim <- c(0,2600) ##Peak
      swe.axis <- seq(0,2500,500) ##Peak
      ylim <- c(0,2000) ##Winter Mean
      swe.axis <- seq(0,2000,500) ##Mean
      ylim <- c(0,2600) ##Spring Mean
      swe.axis <- seq(0,2500,500) ##Mean
   }
   if (j %in% c(3,4)) {
      ylim <- c(0,2000) ##Peak
      swe.axis <- seq(0,1500,500) ##Peak
      ylim <- c(0,1300) ##Mean
      swe.axis <- seq(0,1000,250) ##Mean
      ylim <- c(0,2000) ##Spring Mean
      swe.axis <- seq(0,1750,250) ##Mean
   }
   if (j %in% c(5,6)) {
      ylim <- c(0,1100)
      swe.axis <- seq(0,800,200)
      ylim <- c(0,1500) ##Spring Mean
      swe.axis <- seq(0,1250,250)
   }
   if (j %in% c(7,8)) {
      ylim <- c(0,700) 
      swe.axis <- seq(0,600,200)
      ylim <- c(0,1000) ##Spring
      swe.axis <- seq(0,800,200)
   }
browser()
   make_swe_vs_temp_anomaly_plot(swe.series,tas.anoms,
                                 tas.axis=seq(-4,8,2),swe.axis=swe.axis,
                                 xlim=c(-5,9),ylim=ylim,
                                 whistler.elevs[j],x.axis=x.axis,y.axis=y.axis)

}

##mtext("Annual Peak SWE (mm)",side=2,outer=TRUE,cex=2.0,line=5.0)
##mtext("Annual Temperature Anomaly (\u00B0C)",side=1,outer=TRUE,cex=2.0,line=5.0)

mtext("Spring Average SWE (mm)",side=2,outer=TRUE,cex=2.0,line=5.0)
mtext("Spring Temperature Anomaly (\u00B0C)",side=1,outer=TRUE,cex=2.0,line=5.0)

dev.off()

##---------------------------------------------------------------------
