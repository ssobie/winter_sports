##Script to plot time series of frost free days for Vancouver Intl.

library(ncdf4)
library(PCICt)
library(rgdal)
library(rgeos)
library(zoo)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/data/projects/rci/bcgov/moti/nrcan-precip_case_studies/code/moti.climdex.robjects.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/station.snow.ground.r') 

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
  print(lon.ix)
  print(lat.ix)
  
  data <- matrix(NA,nrow=151,ncol=length(gcm.list))
  yrs <- 1950:2100

  if (grepl('ann',type)) {
      flag <- '*'
  }
  if (grepl('seas',type)) {
      flag <- switch(ix,DJF='-01-',
                        MAM='-04-',
                        JJA='-07-',
                        SON='-10-')
  }
  if (grepl('mon',type)) {
      flag <- switch(ix,DJF='(-12-|-01-|-02-)',
                        MAM='(-03-|-04-|-05-)',
                        JJA='(-06-|-07-|-08-)',
                        SON='(-09-|-10-|-12-)')
  }
         
  for (g in seq_along(gcm.list)) {
    gcm <- gcm.list[g]
    print(gcm)
    var.file <- list.files(path=paste0(read.dir,gcm),pattern=paste0(var.name,'_',type),full.name=TRUE)
    nc <- nc_open(var.file)
    var.series <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    var.time <- netcdf.calendar(nc)   

    ###----------------------------------------------
    ##Snow free days
    snow.yr.fac <- as.factor(format(var.time,'%Y'))
    snow.mn.fac <- as.factor(format(var.time,'%m'))
    seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
    seasonal.fac <- factor(seasons[snow.mn.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))
    avg.fac <- list(snow.yr.fac,seasonal.fac)
    snow.days <- (var.series*100) < 5
    snow.seas.sums <- tapply(snow.days,avg.fac,sum,na.rm=T)
    series.ix <- yrs %in% as.numeric(levels(snow.yr.fac))
    data[series.ix,g] <- snow.seas.sums[,1]
    ###----------------------------------------------
    if (1==0) {

    time.ix <- grep(flag,var.time)
 
    var.sub <- var.series[time.ix]
    var.yrs <- as.numeric(format(var.time,'%Y'))
    series.ix <- yrs %in% var.yrs[time.ix] 
    data[series.ix,g] <- var.sub
    }
    nc_close(nc)
  }
  return(data)
}

##---------------------------------------------------------------------

rcp85.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G','HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')
seas <- 'DJF'
lonc <- -122.95
latc <- 50.13

var.name <- 'snowdepth'

###rcp85.data <- read.gcm.cell(rcp85.list,var.name,'seasonal_mean_BCCAQ',seas,'rcp85',lonc,latc)
rcp85.data <- read.gcm.cell(rcp85.list,var.name,'BCCAQ2',seas,'rcp85',lonc,latc)
var.name <- 'snowdepth'

browser()
##---------------------------------------------------------------------
flag <- switch(seas,DJF='-01-',MAM='-04-',JJA='-07-',SON='-10-',Ann='*')
##flag <- '*'

##ERA-Interim Data
era.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snowdepth_seasonal_mean_BCCAQ2-PRISM_ERA_19790101-20181031.nc'
era.nc <- nc_open(era.file)         
era.lon <- ncvar_get(era.nc,'lon')
era.lat <- ncvar_get(era.nc,'lat')
elon.ix <- which.min(abs(lonc-era.lon))
elat.ix <- which.min(abs(latc-era.lat))
era.series <- ncvar_get(era.nc,var.name,start=c(elon.ix,elat.ix,1),count=c(1,1,-1))
era.time <- netcdf.calendar(era.nc)
nc_close(era.nc)

era.ix <- grep(flag,era.time)
era.sub <- era.series[era.ix][-1]


##NCEP2 Data
ncep2.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snowdepth_seasonal_mean_BCCAQ2-PRISM_NCEP2_19790101-20181031.nc'
ncep2.nc <- nc_open(ncep2.file)        
ncep2.lon <- ncvar_get(ncep2.nc,'lon')
ncep2.lat <- ncvar_get(ncep2.nc,'lat')
nlon.ix <- which.min(abs(lonc-ncep2.lon))
nlat.ix <- which.min(abs(latc-ncep2.lat))
ncep2.series <- ncvar_get(ncep2.nc,var.name,start=c(nlon.ix,nlat.ix,1),count=c(1,1,-1))
 
ncep2.time <- netcdf.calendar(ncep2.nc)
nc_close(ncep2.nc)

ncep2.ix <- grep(flag,ncep2.time)
ncep2.sub <- ncep2.series[ncep2.ix][-1]

##---------------------------------------------------------------------
##Snow Pillow Data
pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/tenquille_lake_asp.csv',sep='')
pillow.data <- read.csv(pillow.file,header=T,as.is=T)
pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
pillow.years <- format(as.Date(pillow.data[,2]),'%Y')
pillow.swe <- pillow.data[,11] ##mm
pillow.peak.swe <- tapply(pillow.swe,as.factor(pillow.years),max,na.rm=T)

##---------------------------------------------------------------------

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

##*********************************************************************

##save(rcp85.data,file='/storage/data/projects/rci/data/winter_sports/ncc_2019/data_files/MV_watersheds.gcms.RData')


##load('/storage/data/projects/rci/data/winter_sports/ncc_2019/data_files/MV_watersheds.gcms.RData')

rcp85.series <- apply(rcp85.data,1,mean,na.rm=T)

ix <- 52:140
px <- 1:52
rx <- 11
yrs <- 2007:2095
hys <- 1956:2007

rcp85.mean <- rollmean(rcp85.series,rx)

hist.90 <- rollmean(apply(rcp85.data,1,quantile,0.9,na.rm=T),rx)[px]
hist.10 <- rollmean(apply(rcp85.data,1,quantile,0.1,na.rm=T),rx)[px]

rcp85.90 <- rollmean(apply(rcp85.data,1,quantile,0.9,na.rm=T),rx)[ix]
rcp85.10 <- rollmean(apply(rcp85.data,1,quantile,0.1,na.rm=T),rx)[ix]

plot.dir <- '/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/'
plot.file <- paste0(plot.dir,'whistler.village.series.winter.mean.depth.clims.png')

##png(plot.file,width=1200,height=900)
png(file=plot.file,width=7,height=3,units='in',res=600,pointsize=6,bg='white')
##pdf(file=plot.file,width=9,height=6,bg='white')

par(mar=c(5,5,5,3))
plot(1951:2100,rcp85.series[-1]*100,type='l',lwd=4,col='white',ylim=c(0,200),
main='Whistler Village Winter Mean Snowdepth (670m)',xlab='Year',ylab='Snowdepth (cm)',xaxs='i',yaxs='i',
cex.axis=1.75,cex.lab=1.75,cex.main=2)
##polygon(c(yrs,rev(yrs)),c(rcp85.10,rev(rcp85.90)),col=alpha('blue',0.3),border=alpha('blue',0.2))
##polygon(c(hys,rev(hys)),c(hist.10,rev(hist.90)),col=alpha('gray',0.3),border=alpha('gray',0.5))

apply(rcp85.data[-1,]*100,2,function(y,x){lines(x,y,col='lightblue',lwd=0.5)},1951:2100)
##lines(1980:2018,ncep2.sub*100,lwd=1,col='darkgreen')
lines(1980:2018,era.sub*100,lwd=1,col='darkgreen')
lines(1951:2100,rcp85.series[-1]*100,lwd=1,col='blue')

##lines(1980:2018,snow.peak*4.2,col='red',lwd=1)
lines(1980:2018,snow.djf,col='red',lwd=1)

#lines(as.numeric(unique(snodas.years)),snodas.peak.swe,lwd=1,col='red')
#lines(as.numeric(unique(pillow.years))[-1],pillow.peak.swe[-1],lwd=1,col='black')


##---------------------------
##Climatologies
if (1==1) {
lines(c(1971,2000),rep(mean(rcp85.series[22:51]*100),2),col='blue',lwd=3)
lines(c(2041,2070),rep(mean(rcp85.series[92:121]*100),2),col='blue',lwd=3)
lines(c(2071,2100),rep(mean(rcp85.series[122:151]*100),2),col='blue',lwd=3)

print(mean(rcp85.series[22:51]*100))
print(mean(rcp85.series[92:121]*100))
print(mean(rcp85.series[122:151]*100))

text(1985,150,'67 cm',cex=1.75)
text(2055,100,'21 cm',cex=1.75)
text(2085,70,'8 cm',cex=1.75)
}
##---------------------------

##lines(hys,rcp85.mean[px],lwd=2,col='darkgray')

##abline(h=seq(-2,8,2),col='gray',lty=3,lwd=1)

box(which='plot')

##legend('topright',legend=c('ERAI','RCP8.5'),col=c('darkgreen','blue'),cex=1.15,pch=15)
legend('topright',legend=c('Station','ERAI','RCP8.5'),col=c('red','darkgreen','blue'),cex=1.15,pch=15)
dev.off()

browser()

rcp26.mean <- rollmean(rcp26.series,11)
rcp45.mean <- rollmean(rcp45.series,11)
rcp85.mean <- rollmean(rcp85.series,11)


plot.file <- paste0(plot.dir,'nfld.annual.tas.smoothed.png')
##png(plot.file,width=900,height=900)
par(mar=c(5,5,5,3))
plot(1955:2095,rcp26.mean,type='l',lwd=4,col='green',ylim=c(0,15),
     main='NFLD Smoothed Annual Average Temperatures',xlab='Year',ylab='TAS (degC)',
     cex.axis=2,cex.lab=2,cex.main=2.5)
lines(1955:2095,rcp45.mean,lwd=4,col='orange')
lines(1955:2094,rcp85.mean,lwd=4,col='red')
abline(h=seq(0,20,5),col='gray',lty=3,lwd=3)
legend('topleft',legend=c('RCP8.5','RCP4.5','RCP2.6'),col=c('red','orange','green'),cex=2,pch=15)
box(which='plot')
##dev.off()
