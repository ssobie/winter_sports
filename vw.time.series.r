##Script to plot time series of frost free days for Vancouver Intl.

library(ncdf4)
library(PCICt)
library(rgdal)
library(rgeos)
library(zoo)
library(scales)
library(raster)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/data/projects/rci/bcgov/moti/nrcan-precip_case_studies/code/moti.climdex.robjects.r',chdir=T)

read.gcm.data <- function(gcm.list,var.name,type,ix,scenario,clip.shp) {

  read.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/')
  
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
    var.raw <- gcm.netcdf.climatologies(var.file,var.name,gcm,NULL,clip.shp)
    var.raw[var.raw > 4] <- NA
    var.series <- apply(var.raw,2,mean,na.rm=T)
    nc <- nc_open(var.file)
    var.time <- netcdf.calendar(nc)   
    time.ix <- grep(flag,var.time)
 
    var.sub <- var.series[time.ix]
    var.yrs <- as.numeric(format(var.time,'%Y'))
    series.ix <- yrs %in% var.yrs[time.ix] 

    data[series.ix,g] <- var.sub
    nc_close(nc)
  }
  return(data)
}

##---------------------------------------------------------------------

rcp85.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G','HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')
seas <- 'MAM'

clip.shp <- readOGR('/storage/data/projects/rci/data/assessments/metro_van/shapefiles/','MVWaterSheds', stringsAsFactors=F)
###clip.shp <- readOGR('/storage/data/projects/rci/data/assessments/whistler/shapefiles/','WhistlerLandscapeUnit', stringsAsFactors=F)
var.name <- 'swe'

rcp85.data <- read.gcm.data(rcp85.list,var.name,'seasonal_mean_BCCAQ',seas,'rcp85',clip.shp)

##---------------------------------------------------------------------
flag <- switch(seas,DJF='-01-',MAM='-04-',JJA='-07-',SON='-10-',Ann='*')
##flag <- '*'

##ERA-Interim Data
era.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/swe_seasonal_mean_BCCAQ2-PRISM_ERA_19790101-20181031.nc'
era.nc <- nc_open(era.file)         
era.time <- netcdf.calendar(era.nc)
nc_close(era.nc)

era.ix <- grep(flag,era.time)
era.raw <- gcm.netcdf.climatologies(era.file,var.name,'ERA',NULL,clip.shp)
era.raw[era.raw>4] <- NA
era.series <- apply(era.raw,2,mean,na.rm=T)
era.sub <- era.series[era.ix]


##NCEP2 Data
ncep2.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/swe_seasonal_mean_BCCAQ2-PRISM_NCEP2_19790101-20181031.nc'
ncep2.nc <- nc_open(ncep2.file)         
ncep2.time <- netcdf.calendar(ncep2.nc)
nc_close(ncep2.nc)

ncep2.ix <- grep(flag,ncep2.time)
ncep2.raw <- gcm.netcdf.climatologies(ncep2.file,var.name,'NCEP2',NULL,clip.shp)
ncep2.raw[ncep2.raw > 4] <- NA
ncep2.series <- apply(ncep2.raw,2,mean,na.rm=T)
ncep2.sub <- ncep2.series[ncep2.ix]

##---------------------------------------------------------------------

##SNODAS Cell
snodas.dir <- '/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/'
snodas.file <- 'swe_seasonal_mean_snodas_modis_grid_van_whistler_20100101-20181231.nc'
snc <- nc_open(paste0(snodas.dir,snodas.file))
snodas.time <- netcdf.calendar(snc)
snodas.ix <- grep(flag,snodas.time)
snodas.series <- apply(gcm.netcdf.climatologies(paste0(snodas.dir,snodas.file),var.name,'SNODAS',NULL,clip.shp),2,mean,na.rm=T)
snodas.sub <- snodas.series[snodas.ix]
snodas.date.sub <- snodas.time[snodas.ix]
nc_close(snc)

##---------------------------------------------------------------------

##save(rcp85.data,file='/storage/data/projects/rci/data/winter_sports/ncc_2019/data_files/MV_watersheds.gcms.RData')
##browser()

##load('/storage/data/projects/rci/data/winter_sports/ncc_2019/data_files/MV_watersheds.gcms.RData')

rcp85.series <- apply(rcp85.data,1,mean,na.rm=T)

ix <- 52:140
px <- 1:52
rx <- 11
yrs <- 2007:2095
hys <- 1956:2007

rcp85.mean <- rollmean(rcp85.series,rx)


##Find the frequency with which back-to-back years of
##low snow levels occurs.

low.snow <- matrix(NA,nrow=151,ncol=12)
rep.lows <- matrix(NA,nrow=150,ncol=12)
for (j in 1:12) {
   q10 <- 0.10*mean(rcp85.data[22:51,j]*1000)
   low.snow[,j] <- rcp85.data[,j]*1000 < q10
   rep.lows[,j] <- (low.snow[2:151,j] == TRUE) & (low.snow[1:150,j] == TRUE)
}

plot.dir <- '/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/'

plot.file <- paste0(plot.dir,'mv.watersheds.repeat.10.percent.lows.spring.mean.swe.bars.png')

png(file=plot.file,width=7,height=3,units='in',res=600,pointsize=6,bg='white')

par(mar=c(5,5,5,3))
plot(1951:2100,apply(rep.lows,1,sum)/12,type='h',lwd=3,col='blue',ylim=c(0,1),lend=2,
main='Metro Vancouver Watersheds Repeat Low SWE',xlab='Year',ylab='Fraction',xaxs='i',yaxs='i',
cex.axis=1.75,cex.lab=1.75,cex.main=2)
##polygon(c(yrs,rev(yrs)),c(rcp85.10,rev(rcp85.90)),col=alpha('blue',0.3),border=alpha('blue',0.2))
##polygon(c(hys,rev(hys)),c(hist.10,rev(hist.90)),col=alpha('gray',0.3),border=alpha('gray',0.5))

box(which='plot')
legend('topleft',legend=c('RCP8.5'),col=c('blue'),cex=1.15,pch=15)
dev.off()

browser()

hist.90 <- rollmean(apply(rcp85.data,1,quantile,0.9,na.rm=T),rx)[px]
hist.10 <- rollmean(apply(rcp85.data,1,quantile,0.1,na.rm=T),rx)[px]

rcp85.90 <- rollmean(apply(rcp85.data,1,quantile,0.9,na.rm=T),rx)[ix]
rcp85.10 <- rollmean(apply(rcp85.data,1,quantile,0.1,na.rm=T),rx)[ix]

plot.dir <- '/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/'

plot.file <- paste0(plot.dir,'mv.watersheds.series.peak.swe.clims.png')
clims <- TRUE

##png(plot.file,width=1200,height=900)
png(file=plot.file,width=7,height=3,units='in',res=600,pointsize=6,bg='white')
##pdf(file=plot.file,width=9,height=6,bg='white')

par(mar=c(5,5,5,3))
plot(1951:2100,rcp85.series[-1],type='l',lwd=4,col='white',ylim=c(0,2100),
main='Metro Vancouver Watersheds Peak SWE',xlab='Year',ylab='SWE (mm)',xaxs='i',yaxs='i',
cex.axis=1.75,cex.lab=1.75,cex.main=2)
##polygon(c(yrs,rev(yrs)),c(rcp85.10,rev(rcp85.90)),col=alpha('blue',0.3),border=alpha('blue',0.2))
##polygon(c(hys,rev(hys)),c(hist.10,rev(hist.90)),col=alpha('gray',0.3),border=alpha('gray',0.5))

if (clims) {
   lwd <- 1
}else{
   lwd <- 2
}
print(lwd)
apply(rcp85.data[-1,]*1000,2,function(y,x){lines(x,y,col='lightblue',lwd=lwd/2)},1951:2100)
lines(1980:2018,era.sub[-1]*1000,lwd=lwd,col='darkgreen')
lines(as.numeric(format(snodas.date.sub,'%Y')),snodas.sub,lwd=lwd,col='red')
##lines(1980:2018,ncep2.sub[-1]*1000,lwd=lwd,col='darkgreen')
lines(1951:2100,rcp85.series[-1]*1000,lwd=lwd,col='blue')


##---------------------------
##Climatologies
if (clims) {
lines(c(1971,2000),rep(mean(rcp85.series[22:51]*1000),2),col='blue',lwd=3)
lines(c(2041,2070),rep(mean(rcp85.series[92:121]*1000),2),col='blue',lwd=3)
lines(c(2071,2100),rep(mean(rcp85.series[122:151]*1000),2),col='blue',lwd=3)

print(mean(rcp85.series[22:51]*1000))
print(mean(rcp85.series[92:121]*1000))
print(mean(rcp85.series[122:151]*1000))

text(1985,1250,'690 mm',cex=1.75)
text(2055,700,'245 mm',cex=1.75)
text(2085,300,'110 mm',cex=1.75)
}
##---------------------------



box(which='plot')
legend('topright',legend=c('ERAI','SNODAS','RCP8.5'),col=c('darkgreen','red','blue'),cex=1.15,pch=15)
dev.off()

