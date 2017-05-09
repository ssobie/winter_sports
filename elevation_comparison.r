##Script to plot the fraction of useable MODIS data
library(ncdf4)
library(PCICt)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

model <- 'ERA'

low.ix <- c(130,65)
high.ix <- c(141,106)
coords <- low.ix
##---------------------------------------------------------------------------------
##MODIS 
modis.dir <- '/storage/data/projects/rci/data/winter_sports/MODIS_VAN_WHISTLER/'
snc.file <- paste0(modis.dir,'snc.modis.van_whistler.20010101-20151231.nc')
snc.nc <- nc_open(snc.file)
modis.time <- netcdf.calendar(snc.nc)

##---------------------------------------------------------------------------------
##SNOW MODEL
snow.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/',model,'/')
snw.file <- paste0(snow.dir,'snowdepth_BCCAQ-PRISM_',model,'_rcp85_r1_1979-2016.nc')
snw.nc <- nc_open(snw.file)
snow.time <- netcdf.calendar(snw.nc)
date.match <- format(snow.time,'%Y-%m-%d') %in% format(modis.time,'%Y-%m-%d')
bnds <- range(which(date.match))


snc.modis <- ncvar_get(snc.nc,'snc',start=c(coords,1),count=c(1,1,-1))
snc.modis[snc.modis > 0 & snc.modis <=100] <- 1
snc.modis[snc.modis > 100] <- NA
snc.valid <- snc.modis

snw.data <- ncvar_get(snw.nc,'snowdepth',start=c(coords,bnds[1]),count=c(1,1,diff(bnds)+1))
snw.cover <- snw.data
##snw.cover[snw.cover>0] <- 1
snw.cover[snw.cover<=0.01] <- 0
  
snw.flagged <- snw.cover
snw.flagged[is.na(snc.valid)] <- NA
snw.check <- snw.flagged
snw.check[snw.check>0] <- 1

snow.diff <- snw.check-snc.valid
snow.upper <- snow.lower <- snow.diff
snow.upper[snow.upper < 1] <- NA
snow.lower[snow.lower > -1] <- NA

subset <- 1:1000
plot(snw.cover[subset],type='l',lwd=3,col='blue',ylim=c(-1,1))
points(snc.valid[subset]*0.1,pch=15,col='black')
points(snow.upper[subset],pch=24,col='green')
points(snow.lower[subset],pch=25,col='red')


nc_close(snc.nc)
nc_close(snw.nc)


