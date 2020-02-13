##Script to plot time series of frost free days for Vancouver Intl.

library(ncdf4)
library(PCICt)
library(rgdal)
library(rgeos)
library(zoo)
library(scales)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

get.coordinates <- function(site) {

  coordinates <- list(spuzzum_creek=c(-121.686,49.674,1197),
                      chilliwack_river=c(-121.71667,49.0333,1600),
                      upper_squamish=c(-123.4333,50.1500,1340),
                      tenquille_lake=c(-122.9333,50.5333,1680))
  rv <- coordinates[[site]]
  return(rv)
}

##-----------------------------------------------------------------

get_erai_data <- function(lonc,latc) {

   ##ERA-Interim Data
   era.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/swe_annual_maximum_BCCAQ2-PRISM_ERA_19790101-20181031.nc'
   era.nc <- nc_open(era.file)         
   era.lon <- ncvar_get(era.nc,'lon')
   era.lat <- ncvar_get(era.nc,'lat')
   elon.ix <- which.min(abs(lonc-era.lon))
   elat.ix <- which.min(abs(latc-era.lat))
   era.series <- ncvar_get(era.nc,var.name,start=c(elon.ix,elat.ix,1),count=c(1,1,-1))*1000
   era.time <- format(netcdf.calendar(era.nc),'%Y-%m-%d')
   nc_close(era.nc)

   return(list(time=era.time,series=era.series))   
}

##-----------------------------------------------------------------

get_ncep2_data <- function(lonc,latc) {
   ##NCEP2 Data
   ncep2.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/swe_annual_maximum_BCCAQ2-PRISM_NCEP2_19790101-20181031.nc'
   ncep2.nc <- nc_open(ncep2.file)        
   ncep2.lon <- ncvar_get(ncep2.nc,'lon')
   ncep2.lat <- ncvar_get(ncep2.nc,'lat')
   nlon.ix <- which.min(abs(lonc-ncep2.lon))
   nlat.ix <- which.min(abs(latc-ncep2.lat))
   ncep2.series <- ncvar_get(ncep2.nc,var.name,start=c(nlon.ix,nlat.ix,1),count=c(1,1,-1))*1000
   ncep2.time <- netcdf.calendar(ncep2.nc)
   nc_close(ncep2.nc)
   return(list(time=ncep2.time,series=ncep2.series))   
}

##-----------------------------------------------------------------

get_snow_pillow_obs <- function(site) {

   ##Snow Pillow Data
   pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'_asp.csv',sep='')
   pillow.data <- read.csv(pillow.file,header=T,as.is=T)
   pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
   pillow.years <- format(as.Date(pillow.data[,2]),'%Y')
   pillow.swe <- pillow.data[,11] ##mm
   pillow.peak.swe <- tapply(pillow.swe,as.factor(pillow.years),max,na.rm=T)
   peak.years <- as.Date(paste0(levels(as.factor(pillow.years)),'-01-01'))

   return(list(series=pillow.peak.swe,time=peak.years))

}

get_snodas_obs <- function(lonc,latc) {

   ##SNODAS Cell
   snodas.dir <- '/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/'
   snodas.file <- 'swe_snodas_modis_grid_van_whistler_20100101-20181231.nc'
   snc <- nc_open(paste0(snodas.dir,snodas.file))
   lon <- ncvar_get(snc,'lon')
   lat <- ncvar_get(snc,'lat')
   lon.ix <- which.min(abs(lonc-lon))
   lat.ix <- which.min(abs(latc-lat))

   snodas.dates <- as.character(netcdf.calendar(snc))
   snodas.years <- format(as.Date(snodas.dates),'%Y')
   snodas.swe <- ncvar_get(snc,'swe',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
   nc_close(snc)

   snodas.peak.swe <- tapply(snodas.swe,as.factor(snodas.years),max,na.rm=T)
   peak.years <- as.Date(paste0(levels(as.factor(snodas.years)),'-01-01'))
   return(list(series=snodas.peak.swe,time=peak.years))
}

##-----------------------------------------------------------------

read.gcm.cell <- function(gcm.list,var.name,type,ix,scenario,lonc,latc) {

  read.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/')

  ##Coordinate indices
  coord.file <- list.files(path=paste0(read.dir,'ACCESS1-0'),pattern=paste0(var.name,'_',type),full.name=TRUE)
  cnc <- nc_open(coord.file)
  lon <- ncvar_get(cnc,'lon')
  lat <- ncvar_get(cnc,'lat')
  nc_close(cnc)
  lon.ix <- which.min(abs(lonc-lon))
  lat.ix <- which.min(abs(latc-lat))
  print('Coordinates indices')
  
  data <- matrix(NA,nrow=151,ncol=length(gcm.list))
  yrs <- 1950:2100
  dates <- seq(from=as.Date('1950-01-01'),by='year',to=as.Date('2100-12-01'))
  for (g in seq_along(gcm.list)) {
    gcm <- gcm.list[g]
    print(gcm)
    var.file <- list.files(path=paste0(read.dir,gcm),pattern=paste0(var.name,'_',type),full.name=TRUE)
    nc <- nc_open(var.file)
    var.series <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    var.time <- netcdf.calendar(nc)   
    nc_close(nc)

    var.sub <- var.series
    var.yrs <- as.numeric(format(var.time,'%Y'))
    series.ix <- yrs %in% var.yrs
    data[series.ix,g] <- var.sub*1000

  }
  return(list(series=data,time=dates))
}

##---------------------------------------------------------------------

rcp85.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G','HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')
seas <- 'Annual'
var.name <- 'swe'

sites <- c('spuzzum_creek','upper_squamish','chilliwack_river','tenquille_lake')
site.names <- c('Spuzzum Creek','Upper Squamish','Chilliwack River','Tenquille Lake')

plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
png(file=paste0(plot.dir,'future.SWE.pillow.series.2019.png'),width=6,height=5,units='in',res=600,pointsize=6,bg='white')
par(mfrow=c(4,1),oma=c(1,2,1,1))

##Loop over sites

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)
    coords <- get.coordinates(site)
    lonc <- coords[1]
    latc <- coords[2]
    rcp85.data <- read.gcm.cell(rcp85.list,var.name,'annual_maximum',seas,'rcp85',lonc,latc)
    rcp85.series <- apply(rcp85.data$series,1,mean,na.rm=T)
    ix <- 52:140
    px <- 1:52
    rx <- 11
    yrs <- 2007:2095
    hys <- 1956:2007
    rcp85.mean <- rollmean(rcp85.series,rx)
    
    erai.peak <- get_erai_data(lonc,latc)
    ncep2.peak <- get_ncep2_data(lonc,latc)

    pillow.peak <- get_snow_pillow_obs(site)
    snodas.peak <- get_snodas_obs(lonc,latc)

    par(mar=c(4,6,2,2))
    plot(rcp85.data$time,rcp85.series,xlab='Date',ylab='SWE (mm)',yaxs='i',
           type='l',lwd=2.5,col='blue',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,xaxs='i',
           xlim=c(as.Date('1951-01-01'),as.Date('2100-01-31')),ylim=c(0,3200))
      apply(rcp85.data$series,2,function(x,y){lines(y,x,col='lightblue',lwd=2.0)},rcp85.data$time)
      points(as.Date(snodas.peak$time),snodas.peak$series,col=alpha('red',0.5),lwd=1)
      points(as.Date(pillow.peak$time),pillow.peak$series,pch=16,col='black')
      points(as.Date(erai.peak$time),erai.peak$series,pch=16,col='green')
      lines(rcp85.data$time,rcp85.series,col='blue',lwd=2.5)
      text(as.Date('2080-01-01'),2500,site.names[i],cex=2.5)
      if (i==4) {
         legend('bottomleft',legend=c('ASP Obs.','SNODAS','ERA','NCEP'),
                             col=c('black','red','lightblue','blue'),pch=16,cex=1.75)
      }
    box(which='plot')
}

dev.off()
##---------------------------------------------------------------------



