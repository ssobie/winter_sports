##Script to plot time series of frost free days for Vancouver Intl.

library(ncdf4)
library(PCICt)
library(rgdal)
library(rgeos)
library(zoo)
library(raster)
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

read_cell_clim <- function(gcm,var.name,var.file,read.dir,lon.ix,lat.ix,seas.ix) { ##,lonc,latc) {

  ##Coordinate indices
  nc <- nc_open(paste0(read.dir,var.file))
  var.clim <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,seas.ix),count=c(1,1,1))
  var.time <- format(netcdf.calendar(nc),'%Y-%m-%d')

  nc_close(nc)
  return(var.clim)
}

##---------------------------------------------------------------------

read_cell_series <- function(gcm,var.name,var.file,read.dir,cell.ix,seas) { ##,lonc,latc) {

  ##Coordinate indices
  nc <- nc_open(paste0(read.dir,var.file))
  var.time <- format(netcdf.calendar(nc),'%Y-%m-%d')
  seas.ix <- switch(seas,
                    winter='-01-',
                    spring='-04-',
                    summer='-07-',
                    fall='-10-',
                    annual='*')
  time.ix <- grep(seas.ix,var.time,'%m')
  time.sub <- var.time[time.ix]
  if (grepl('HadGEM',var.file)) {
     years <- 1950:2100
     had.years <- as.numeric(format(as.Date(time.sub),'%Y'))
     had.ix <- years %in% had.years
     had.fill <- rep(NA,length(years))
     had.time <- paste0(years,'-',format(as.Date(time.sub[20]),'%m-%d'))
     time.sub <- had.time
  }

  var.matrix <- matrix(NA,nrow=151,ncol=nrow(cell.ix))

  for (j in 1:nrow(cell.ix)) { 
     print(j)
     var.series <- ncvar_get(nc,var.name,start=c(cell.ix[j,1],cell.ix[j,2],1),count=c(1,1,-1))
     var.sub <- var.series[time.ix]
     if (grepl('HadGEM',var.file)) {
        had.fill[had.ix] <- var.sub
        var.sub <- had.fill
     }
     var.matrix[,j] <- var.sub
  }

  nc_close(nc)
  rv <- var.matrix  ##list(series=var.matrix,time=time.sub)
  return(rv)
}

##---------------------------------------------------------------------
##*********************************************************************

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

prism.file <- '/storage/data/projects/rci/data/prism/van_whistler_prism_dem_elevations.nc'
dem.brick <- brick(prism.file)

##Indices of cells within Whistler-Blackcomb Area

shape.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'
clip.shp <- readOGR(shape.dir,'whistler_blackcomb_buffer')
clip.84 <- spTransform(clip.shp,CRS("+init=epsg:4326"))
dem.wb <- mask(dem.brick,clip.84)
dem.rotated <- (t(as.matrix(subset(dem.wb,1))[323:1,]))
ix <- which(!is.na(dem.rotated),arr.ind=T)

##------------------------------------------------------------
type <- 'swe_seasonal_mean'
seas.ix <- 2 
seas <- 'spring'

##x11()
##plot(c(),xlim=c(-6,10),ylim=c(500,2600),pch=16,col='darkgray',
##        xlab='TAS Anoms',ylab='Elevation (m)',yaxs='i',cex=1.25)

dem.vals <- dem.rotated[ix]

swe.models <- vector(mode='list',length=length(rcp85.list))
tas.models <- vector(mode='list',length=length(rcp85.list))

for (i in seq_along(rcp85.list)) {
   gcm <- rcp85.list[i]
   print(gcm)
   swe.dir <- paste0(base.dir,'snow_model/calibrated_',gcm,'_PNWNAmet_prism_tps/')
   swe.file <- list.files(path=swe.dir,pattern=type)
   swe.data <- read_cell_series(gcm,var.name='swe',swe.file,swe.dir,ix,seas)*1000
   swe.models[[i]] <- swe.data

   tas.dir <- paste0(base.dir,gcm,'/tas_climatology/')
   tasmax.file <- list.files(path=tas.dir,pattern=paste0('tasmax_seasonal_mean_gcm_prism_',gcm))
   tasmax.data <- read_cell_series(gcm,var.name='tasmax',tasmax.file,tas.dir,ix,seas)
   tasmin.file <- list.files(path=tas.dir,pattern=paste0('tasmin_seasonal_mean_gcm_prism_',gcm))
   tasmin.data <- read_cell_series(gcm,var.name='tasmin',tasmin.file,tas.dir,ix,seas)

   tas.data <- (tasmax.data + tasmin.data) / 2
   tas.anoms <- apply(tas.data,2,function(x){x-mean(x[31:61])})   
   tas.models[[i]] <- tas.anoms


}

tas.levels <- dem.vals*0
tas.lengths <- dem.vals*0

for (j in 1:nrow(ix)) {
   swe.matrix <- matrix(0,nrow=151,ncol=12)
   tas.matrix <- matrix(0,nrow=151,ncol=12)    
   for (k in seq_along(rcp85.list)) {
      swe.cell <- swe.models[[k]]
      swe.matrix[,k] <- swe.cell[,j]
      tas.cell <- tas.models[[k]]
      tas.matrix[,k] <- tas.cell[,j]
   }
   swe.thresh <- apply(swe.matrix,2,function(x){return(x > 110 & x < 130)})
   swe.thresh[is.na(swe.thresh)] <- FALSE
   swe.crit <- swe.matrix[swe.thresh]
   tas.crit <- tas.matrix[swe.thresh]
   print('Number of critical values')
   print(length(tas.crit))
   tas.lengths[j] <- length(tas.crit)
   tas.levels[j] <- mean(tas.crit)
   if (quantile(swe.matrix,0.9,na.rm=T,names=F) > 5000) {
      tas.levels[j] <- NA
   }
}

##Use values with at least 3 critical tas values
tx <- tas.lengths > 3

##Make plot of temperature dependence of critical swe levels

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/whistler.min.ski.elevation.vs.tas.png')

png(file=plot.file,width=4,height=4,units='in',res=600,pointsize=6,bg='white')

xlim=c(0,7)
ylim=c(500,2000)

plot(tas.levels[tx],dem.vals[tx],xlim=xlim,ylim=ylim,pch=16,col='darkgray',
        xlab='Spring Temperature Anomaly (\u00B0C)',
        ylab='Min. SWE Elevation (m)',xaxs='i',yaxs='i',cex=1.35,
        cex.axis=1.35,cex.lab=1.35)
   abline(lm(dem.vals[tx]~tas.levels[tx]))
box(which='plot')
dev.off()

##---------------------------------------------------------------------

