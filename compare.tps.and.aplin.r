##Script to compare the differences between ANUSPLIN and TPS
##along with the BCCAQ2 bias corrected versions and swe outputs.

library(ncdf4)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R')

##Coords for the cell north of Spuzzum Creek

lon.c <-  -121.7
lat.c <- 49.75214

compare.series <- function(varname,a.file,t.file,lon.c,lat.c) {
 
   anc <- nc_open(a.file)
   alon <- ncvar_get(anc,'lon')
   alat <- ncvar_get(anc,'lat')
   aix <- c(which.min(abs(lon.c-alon)),which.min(abs(lat.c-alat)))
   print(aix)
   if (a.file=='/storage/data/climate/downscale/BCCAQ2/ANUSPLIN/anusplin_pr_final.nc') {
      load('/storage/data/projects/rci/data/winter_sports/BCCAQ2/ANUSPLIN300/pr_232_106_canada.RData')
   } else {
      apr <- ncvar_get(anc,varname,start=c(aix,1),count=c(1,1,-1))
   }
##save(apr,file='/storage/data/projects/rci/data/winter_sports/BCCAQ2/ANUSPLIN300/pr_232_106_canada.RData') 

   atime <- netcdf.calendar(anc)
   amon.fac <- as.factor(format(atime,'%Y-%m'))
   apr.mon <- tapply(apr,amon.fac,sum)
   amons <- as.Date(paste0(levels(amon.fac),'-01'))

   tnc <- nc_open(t.file)
   tlon <- ncvar_get(tnc,'lon')
   tlat <- ncvar_get(tnc,'lat')
   tix <- c(which.min(abs(lon.c-tlon)),which.min(abs(lat.c-tlat)))
   print(tix)
   tpr <- ncvar_get(tnc,varname,start=c(tix,1),count=c(1,1,-1))
   ttime <- netcdf.calendar(tnc)
   tmon.fac <- as.factor(format(ttime,'%Y-%m'))
   tpr.mon <- tapply(tpr,tmon.fac,sum)
   tmons <- as.Date(paste0(levels(tmon.fac),'-01'))

   ##rv <- list(amons=amons,apr=apr.mon,tmons=tmons,tpr=tpr.mon)
   rv <- list(amons=atime,apr=apr,tmons=ttime,tpr=tpr)
   nc_close(anc)
   nc_close(tnc)
   return(rv)
}

##Observations
a.file <- '/storage/data/climate/downscale/BCCAQ2/ANUSPLIN/anusplin_pr_final.nc'
t.file <- '/storage/data/projects/rci/data/winter_sports/TPS/van.whistler.tps.pr.nc'

obs <- compare.series('pr',a.file,t.file,lon.c,lat.c)

##--
#Compare BCCAQ2 versions of ERA -> ANUSPLIN or TPS

a.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/ERA/pr_day_QDM_ERA_1979-2016.nc'
t.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/ERA/pr_day_QDM_ERA_19790101-20181031.nc'

bccaq2 <- compare.series('pr',a.file,t.file,lon.c,lat.c)


#Compare GCM-PRISM versions of ERA -> BCCAQ2

a.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/ERA/pr_gcm_prism_ERA_1979-2016.nc'
t.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/ERA/pr_gcm_prism_ERA_19790101-20181031.nc'

prism <- compare.series('pr',a.file,t.file,lon.c,lat.c)


par(mfrow=c(3,1))

plot(obs$tmons,obs$tpr,type='l',col='blue',lwd=2,main='Monthly Total Precipitation',ylab='mm',xlab='Date',ylim=c(0,250))
lines(obs$amons,obs$apr,col='red',lwd=2)
abline(h=c(0,50,100,150,200),col='gray',lty=2)
legend('topleft',leg=c('TPS','APLIN'),col=c('blue','red'),pch=15)

plot(bccaq2$tmons,bccaq2$tpr,type='l',col='blue',lwd=2,main='Monthly Total Precipitation',ylab='mm',xlab='Date',ylim=c(0,250))
lines(bccaq2$amons,bccaq2$apr,col='red',lwd=2)
abline(h=c(0,50,100,150,200),col='gray',lty=2)
legend('topleft',leg=c('TPS','APLIN'),col=c('blue','red'),pch=15)

plot(prism$tmons,prism$tpr,type='l',col='white',lwd=2,main='Monthly Total Precipitation',ylab='mm',xlab='Date',ylim=c(0,250))
lines(prism$amons,prism$apr,col='red',lwd=2)
lines(prism$tmons,prism$tpr,col='blue',lwd=2)
abline(h=c(0,50,100,150,200),col='gray',lty=2)
legend('topleft',leg=c('TPS','APLIN'),col=c('blue','red'),pch=15)


##Observations
a.file <- '/storage/data/climate/downscale/BCCAQ2/ANUSPLIN/anusplin_tasmin_final.nc'
t.file <- '/storage/data/projects/rci/data/winter_sports/TPS/van.whistler.tps.tasmin.nc'

obs <- compare.series('tasmin',a.file,t.file,lon.c,lat.c)
plot(obs$tmons,obs$tpr,type='l',col='white',lwd=2,main='Min Temperature ',ylab='mm',xlab='Date')
lines(obs$amons,obs$apr,col='red',lwd=2)
lines(obs$tmons,obs$tpr,col='blue',lwd=2)
abline(h=c(0,50,100,150,200),col='gray',lty=2)
legend('topleft',leg=c('TPS','APLIN'),col=c('blue','red'),pch=15)


##--
##Compare snow model outputs of SWE



