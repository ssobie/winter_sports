##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)

library(ncdf4)
library(PCICt)
library(plotrix)

sites <- c('shovelnose_mountain',
           'brookmere',
           'lightning_lake',
           'callaghan',
           'orchid_lake',
           'palisade_lake',
           'grouse_mountain',
           'dog_mountain',
           'stave_lake',
           'nahatlatch',
           'wahleach',
           'klesilkwa',
           'hamilton_hill',
           'upper_squamish',
           'spuzzum_creek',
           'chilliwack_river',
           'tenquille_lake')

model <- 'ERA5'
site <- 'chilliwack_river'
site.name <- 'Tenquille Lake'

##---------------------------------------------------------------------------------
##Snow Pillow Data
pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'_asp.csv',sep='')
pillow.data <- read.csv(pillow.file,header=T,as.is=T)
pillow.dates <- as.Date(format(as.Date(pillow.data[,2]),'%Y-%m-%d'))
pillow.swe <- pillow.data[,11] ##mm
pillow.pack <- pillow.data[,13] ##cm
pillow.cover <- pillow.swe
pillow.cover[pillow.swe>10] <- 1
pillow.cover[pillow.swe<=10] <- 0


##MODIS 
modis.dir <- '/storage/data/projects/rci/data/winter_sports/MODIS_VAN_WHISTLER/'
mnc.file <- paste0(modis.dir,'snc.modis.van_whistler.20010101-20181231.nc')
mnc.nc <- nc_open(mnc.file)
m.lon <- ncvar_get(mnc.nc,'lon')
m.lat <- ncvar_get(mnc.nc,'lat')
modis.time <- netcdf.calendar(mnc.nc)
modis.dates <- as.Date(format(netcdf.calendar(mnc.nc),'%Y-%m-%d'))

##---------------------------------------------------------------------------------
##SNODAS 
snodas.dir <- '/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/'
snc.file <- paste0(snodas.dir,'swe_snodas_modis_grid_van_whistler_20100101-20181231.nc')
snc.nc <- nc_open(snc.file)
s.lon <- ncvar_get(snc.nc,'lon')
s.lat <- ncvar_get(snc.nc,'lat')
snodas.time <- netcdf.calendar(snc.nc)
snodas.dates <- as.Date(format(netcdf.calendar(snc.nc),'%Y-%m-%d'))
##---------------------------------------------------------------------------------
##SNOW MODEL from Simulations
snow.time.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/calibrated_ERA5_prism_tps/')
esnw.time.file <- paste0(snow.time.dir,'snowdepth_BCCAQ2-PRISM_ERA5_1980-2018.nc')
esnw.nc <- nc_open(esnw.time.file)
era.time <- netcdf.calendar(esnw.nc)
era.dates <- as.Date(format(netcdf.calendar(esnw.nc),'%Y-%m-%d'))

##---------------------------------------------------------------------------------

date.match <- format(era.time,'%Y-%m-%d') %in% format(modis.time,'%Y-%m-%d')
modis.match <- format(modis.time,'%Y-%m-%d') %in% format(era.time,'%Y-%m-%d')
snodas.match <- format(modis.time,'%Y-%m-%d') %in% format(snodas.time,'%Y-%m-%d')

era.snodas.match <- format(era.time,'%Y-%m-%d') %in% format(snodas.time,'%Y-%m-%d')
snodas.era.match <- format(snodas.time,'%Y-%m-%d') %in% format(era.time,'%Y-%m-%d')

modis.era.snodas.match <- format(modis.time,'%Y-%m-%d') %in% format(era.time[era.snodas.match],'%Y-%m-%d')

modis.pillow.match  <- format(modis.time,'%Y-%m-%d') %in% format(pillow.dates,'%Y-%m-%d')
pillow.modis.match  <- format(pillow.dates,'%Y-%m-%d') %in% format(modis.time,'%Y-%m-%d')

bnds <- range(which(date.match))

  snow.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/')
  era.file <- paste0(snow.dir,site,'_ERA5_PRISM_TPS_snow_model_data.csv')
  era.sims <- read.csv(era.file,header=TRUE,as.is=TRUE)
  era.snow <- era.sims$Snowdepth
  era.swe <- era.sims$SWE*1000

  print(site)
  coords <- get_coordinates(site)
  lon <- coords[1]
  elev <- coords[3]
  
  m.ox <- which.min(abs(coords[1]-m.lon))
  m.ax <- which.min(abs(coords[2]-m.lat))

  mnc.data <- ncvar_get(mnc.nc,'snc',start=c(m.ox,m.ax,1),count=c(1,1,-1))[modis.match]
  mnc.modis <- mnc.data
  mnc.modis[mnc.modis > 0 & mnc.modis <=100] <- 1
  mnc.modis[mnc.modis > 100] <- NA
  mnc.valid <- mnc.modis
  ##mnc.modis[mnc.modis == 0]  <- NA

  s.ox <- which.min(abs(coords[1]-s.lon))
  s.ax <- which.min(abs(coords[2]-s.lat))

  snc.data <- ncvar_get(snc.nc,'swe',start=c(s.ox,s.ax,1),count=c(1,1,-1))
  snc.snodas <- snc.data
  snc.snodas[snc.snodas <= 10] <- 0
  snc.snodas[snc.snodas > 10] <- 1

  era.cover <- era.snow ##[bnds[1]:bnds[2],]
  era.cover[era.cover>0.1] <- 1
  era.cover[era.cover<=0.1] <- 0

  ##For snow only

  flag <- is.na(mnc.modis)
  era.snow <- era.cover[bnds[1]:bnds[2]]
  era.snow[flag] <- NA

  mflag <- is.na(mnc.modis[snodas.match])
  modis.snodas <- mnc.modis[snodas.match]
  modis.era.snodas <- mnc.modis[modis.era.snodas.match]
  sflag <- is.na(modis.era.snodas)

  era.snodas <- era.cover[era.snodas.match]
  era.snodas[sflag] <- NA

nc_close(mnc.nc)
nc_close(snc.nc)
nc_close(esnw.nc)


##Create a 3-figure plot
##Time series with snow cover days indicated
##Histogram of Model Days with snow and MODIS without
##Histogram of Model Days without snow and MODIS with snow

layout(mat = rbind(c(1,1),c(2,3)))

      par(mar=c(4,6,2,2))
      plot(era.dates,era.swe,xlab='Date',ylab='SWE (mm)',yaxs='i',
           type='l',lwd=3,col='white',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,xaxs='i',
           xlim=c(as.Date('2010-08-01'),as.Date('2018-12-31')),ylim=c(0,3000))

      ##points(as.Date(pillow.dates),pillow.swe,pch=16,col='black')
      #lines(era.dates,era.swe,col='blue',lwd=2.75)
      text(as.Date('1996-01-01'),2500,site.name,cex=2.5)
      sc.ix <- mnc.valid==1
      abline(v=modis.dates[sc.ix],col='blue',lwd=0.5)
      nc.ix <- mnc.valid==0
      abline(v=modis.dates[nc.ix],col='red',lwd=0.5)
      ##lines(era.dates,era.swe,col='black',lwd=2.75)
      ##points(snodas.dates,snc.data,col='green',lwd=2)
      points(pillow.dates,pillow.swe,col='black',pch=16)
      ##legend('bottomleft',legend=c('ASP Obs.','SNODAS',model),col=c('black','red','blue'),pch=16,cex=1.75)

      me.diff <- (pillow.cover[pillow.modis.match]-mnc.valid[modis.pillow.match])*0.5
      ix <- which(me.diff==0.5)
      ixm <- which(me.diff==-0.5)
      hist(as.numeric(format(modis.dates[ix],'%j')),breaks=seq(1,365,5),col='red',ylim=c(0,20))
      hist(as.numeric(format(modis.dates[ixm],'%j')),breaks=seq(1,365,5),col='blue',ylim=c(0,20))