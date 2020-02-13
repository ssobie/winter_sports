##Script to pull out the time series of pr, tasmax, tamsin from the driving models 
library(ncdf4)
library(raster)


get.coordinates <- function(site) {

  coordinates <- list(callaghan=c(-123.1036,50.1383278,1009),
                      orchid_lake=c(-123.0519638,49.53678,1178),
                      palisade_lake=c(-123.0321944,49.454433,898),
                      grouse_mountain=c(-123.0774472,49.383655,1126),
                      dog_mountain=c(-122.96255,49.37251944,1007),
                      dickson_lake=c(-122.06984166,49.3168194,1147),
                      stave_lake=c(-122.315805,49.58030277,1211),
                      nahatlatch=c(-122.059261,49.825866,1530),
                      wahleach=c(-121.57945,49.2298694,1395),
                      klesilkwa=c(-121.3086527,49.129438,610),
                      lightning_lake=c(-120.850205,49.044788,1254),
                      brookmere=c(-120.87397,49.815027,994),
                      shovelnose_mountain=c(-120.864175,49.8546305,1456),
                      hamilton_hill=c(-120.7955805,49.4988027,1477),
                      spuzzum_creek=c(-121.686,49.74,1197),
                      chilliwack_river=c(-121.71667,49.0333,1600),
                      upper_squamish=c(-123.4333,50.1500,1340),
                      wahleach_lake=c(-121.5833,49.2333,1400),
                      tenquille_lake=c(-122.9333,50.5333,1680))

  rv <- coordinates[[site]]
  return(rv)
}

get_PRISM_values <- function(site, prism.var, file,buffer) {

  coords <- get.coordinates(site)
  plot.title <- site

  nc <- nc_open(file)

  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')

  lon.bnds <- coords[1]
  lat.bnds <- coords[2]

  lon.ix <- which.min(abs(lon-lon.bnds))
  lat.ix <- which.min(abs(lat-lat.bnds))

  prism.values <- ncvar_get(nc,prism.var,start=c(lon.ix-buffer,lat.ix-buffer,1),
                                          count=c(2*buffer+1,2*buffer+1,-1))
  nc_close(nc)
  return(prism.values)
}

get_PRISM_elevs <- function(site,file,buffer) {

  coords <- get.coordinates(site)
  plot.title <- site

  nc <- nc_open(file)

  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')

  lon.bnds <- coords[1]
  lat.bnds <- coords[2]

  lon.ix <- which.min(abs(lon-lon.bnds))
  lat.ix <- which.min(abs(lat-lat.bnds))

  elevations <- ncvar_get(nc,'elevation',start=c(lon.ix-buffer,lat.ix-buffer,1),
                                         count=c(2*buffer+1,2*buffer+1,-1))
  nc_close(nc)
  return(elevations)
}


get_PRISM_tiff_elevs <- function(site, dem, buffer) {

  coords <- get.coordinates(site)
  plot.title <- site

  dem.coords <- coordinates(dem)

  lon <- unique(dem.coords[,1])
  lat <- unique(dem.coords[,2])

  lon.bnds <- coords[1]
  lat.bnds <- coords[2]
  elev <- coords[3]

  lon.ix <- which.min(abs(lon-lon.bnds))
  lat.ix <- which.min(abs(lat-lat.bnds))

  dem.elev <- dem[(lat.ix-buffer):(lat.ix+buffer),(lon.ix-buffer):(lon.ix+buffer),1]
browser()
  return(list(site=elev,dem=dem.elev))
}

  sites <- c('callaghan',
             'orchid_lake',
             'palisade_lake',
             'grouse_mountain',
             'dog_mountain',
             'stave_lake',
             'nahatlatch',
             'wahleach',
             'klesilkwa',
             'lightning_lake',
             'brookmere',
             'shovelnose_mountain',
             'hamilton_hill',
             'spuzzum_creek',
             'chilliwack_river',
             'upper_squamish',
             'tenquille_lake')

##site <- 'palisade_lake'
##site.name <- 'Palisade Lake'
tasmax.lapses <- matrix(0,nrow=length(sites),ncol=12)
tasmin.lapses <- matrix(0,nrow=length(sites),ncol=12)
pr.lapses <- matrix(0,nrow=length(sites),ncol=12)

for (s in seq_along(sites)) {
   site <- sites[s]
   buffer <- 2
   base.dir <- "/storage/data/projects/PRISM/bc_climate/grids/"

##dem.tas <- brick(paste0(base.dir,"PRISM_BC_Domain_30s_dem.grass.tiff"))
##dem.pr <-  brick(paste0(base.dir,"PRISM_BC_Domain_30s375m_dem.grass.test.tiff"))

   prism.dem <- '/storage/data/projects/rci/data/prism/bc_prism_dem_elevations.nc'
   tas.elevs <- get_PRISM_elevs(site, prism.dem,buffer)

   ##PRISM Climatologies
   tasmax.clim.file <- '/storage/data/climate/PRISM/dataportal/tmax_monClim_PRISM_historical_run1_198101-201012.nc'
   tasmin.clim.file <- '/storage/data/climate/PRISM/dataportal/tmin_monClim_PRISM_historical_run1_198101-201012.nc'
   pr.clim.file <- '/storage/data/climate/PRISM/dataportal/pr_monClim_PRISM_historical_run1_198101-201012.nc'

   tasmax.clims <- get_PRISM_values(site,'tmax',tasmax.clim.file,buffer)
   tasmin.clims <- get_PRISM_values(site,'tmin',tasmin.clim.file,buffer)
   pr.clims <- get_PRISM_values(site,'pr',pr.clim.file,buffer)

   elev.order <- tas.elevs[order(as.vector(tas.elevs))]

   for (m in 1:12) {
      tasmax.order <- tasmax.clims[,,m][order(as.vector(tas.elevs))]
      tasmin.order <- tasmin.clims[,,m][order(as.vector(tas.elevs))]
      pr.order <- pr.clims[,,m][order(as.vector(tas.elevs))]
      tasmax.lapses[s,m] <- median(diff(tasmax.order)/diff(elev.order)*100,na.rm=T)
      tasmin.lapses[s,m] <- median(diff(tasmin.order)/diff(elev.order)*100,na.rm=T)
      pr.lapses[s,m] <- median(diff(pr.order)/diff(elev.order)*100,na.rm=T)
   }
}

tasmax.lapses <- rbind(c('Site',month.abb),cbind(sites,round(tasmax.lapses,3)))
tasmin.lapses <- rbind(c('Site',month.abb),cbind(sites,round(tasmin.lapses,3)))
pr.lapses <- rbind(c('Site',month.abb),cbind(sites,round(pr.lapses,3)))

tasmax.file <- '/storage/data/projects/rci/data/winter_sports/course_and_pillow_prism_tasmax_lapse_rates.csv'
tasmin.file <- '/storage/data/projects/rci/data/winter_sports/course_and_pillow_prism_tasmin_lapse_rates.csv'
pr.file <- '/storage/data/projects/rci/data/winter_sports/course_and_pillow_prism_pr_lapse_rates.csv'

write.table(tasmax.lapses,file=tasmax.file,sep=',',col.name=F,row.name=F,quote=F)
write.table(tasmin.lapses,file=tasmin.file,sep=',',col.name=F,row.name=F,quote=F)
write.table(pr.lapses,file=pr.file,sep=',',col.name=F,row.name=F,quote=F)


browser()






plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'

png(file=paste0(plot.dir,site,'.tasmax.vs.elevation.png'),width=5,height=6,units='in',res=300,pointsize=6,bg='white')

par(mfrow=c(4,3))
par(mar=c(4,4,1,1),oma=c(1,1,3,1))
##par(mar=c(0,0,0,0),oma=c(7,7,4,4))
for (m in 1:12) {
   plot(tas.elevs,tasmax.clims[,,m],pch=16,xlab='Elevation (m)',ylab='Temperature (degC)',
        cex=1.5,cex.lab=1.5,cex.axis=1.5)
   if (grepl(m,'(1|4|7|10)')) {
  ##    axis(2,at=0:3,label=0:3)
   }
   if (grepl(m,'(10|11|12)')) {
  ##    axis(1,at=seq(200,1200,200),label=seq(200,1200,200),cex=1.5)
   }
   text(1000,0.9*(par('usr')[4]-par('usr')[3])+par('usr')[3],month.abb[m],cex=2.2)
   mtext(paste0('Maximum Temperature at ',site.name),side=3,cex=1.5,outer=T,line=1)
   box(which='plot')
}

dev.off()


png(file=paste0(plot.dir,site,'.tasmin.vs.elevation.png'),width=5,height=6,units='in',res=300,pointsize=6,bg='white')

par(mfrow=c(4,3))
par(mar=c(4,4,1,1),oma=c(1,1,3,1))
##par(mar=c(0,0,0,0),oma=c(7,7,4,4))
for (m in 1:12) {
   plot(tas.elevs,tasmin.clims[,,m],pch=16,xlab='Elevation (m)',ylab='Temperature (degC)',
        cex=1.5,cex.lab=1.5,cex.axis=1.5)
   if (grepl(m,'(1|4|7|10)')) {
  ##    axis(2,at=0:3,label=0:3)
   }
   if (grepl(m,'(10|11|12)')) {
  ##    axis(1,at=seq(200,1200,200),label=seq(200,1200,200),cex=1.5)
   }
   text(1000,0.9*(par('usr')[4]-par('usr')[3])+par('usr')[3],month.abb[m],cex=2.2)
   mtext(paste0('Minimum Temperature at ',site.name),side=3,line=1,outer=T,cex=1.5)
   box(which='plot')
}

dev.off()

elev.order <- tas.elevs[order(as.vector(tas.elevs))]


png(file=paste0(plot.dir,site,'.tasmax.lapse.rates.png'),width=5,height=6,units='in',res=300,pointsize=6,bg='white')

par(mfrow=c(6,2))
par(mar=c(4,4,1,1),oma=c(1,1,3,1))
##par(mar=c(0,0,0,0),oma=c(7,7,4,4))
for (m in 1:12) {
   tasmax.order <- tasmax.clims[,,m][order(as.vector(tas.elevs))]
   plot(diff(tasmax.order)/diff(elev.order)*100,ylim=c(-15,3),
        pch=18,xlab=month.abb[m],ylab='Lapse Rates (degC / 100m)',
        cex=2.0,cex.lab=1.75,cex.axis=1.75)
   print(median(diff(tasmax.order)/diff(elev.order)*100))
   abline(h=0)
   mtext(paste0('Tasmax Lapse Rates (degC / 100m) at ',site.name),side=3,line=1,outer=T,cex=1.5)
   box(which='plot')
}

dev.off()


png(file=paste0(plot.dir,site,'.tasmin.lapse.rates.png'),width=5,height=6,units='in',res=300,pointsize=6,bg='white')

par(mfrow=c(6,2))
par(mar=c(4,4,1,1),oma=c(1,1,3,1))
##par(mar=c(0,0,0,0),oma=c(7,7,4,4))
for (m in 1:12) {
   tasmin.order <- tasmin.clims[,,m][order(as.vector(tas.elevs))]
   plot(diff(tasmin.order)/diff(elev.order)*100,ylim=c(-15,3),
        pch=18,xlab=month.abb[m],ylab='Lapse Rates (degC / 100m)',
        cex=2.0,cex.lab=1.75,cex.axis=1.75)
   print(median(diff(tasmin.order)/diff(elev.order)*100))
   abline(h=0)
   mtext(paste0('Tasmin Lapse Rates (degC / 100m) at ',site.name),side=3,line=1,outer=T,cex=1.5)
   box(which='plot')
}

dev.off()
