##Script to pull out the time series of pr, tasmax, tamsin from the driving models 
library(ncdf4)
library(raster)

source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)

get_PRISM_values <- function(site, prism.var, file,buffer) {

  coords <- get_coordinates(site)
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

  coords <- get_coordinates(site)
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

  return(list(site=elev,dem=dem.elev))
}


type <- 'evaluation'
sites <- c('blackwall_peak_course','boston_bar_lower','boston_bar_upper','burwell_lake',
           'chapman_creek','cornwall_hills','diamond_head','edwards_lake',
           'hollyburn','hope',
           'loch_lomond','lytton','mount_seymour','new_tashme',
           'ottomite','pavilion_mountain',
           'sumallo_river','tenquille_course','whistler_mountain','wolverine_creek' )

type <- 'calibration'
sites <- c('shovelnose_mountain','brookmere','lightning_lake','callaghan','orchid_lake',
           'palisade_lake','grouse_mountain','dog_mountain','stave_lake','nahatlatch',
           'wahleach','klesilkwa','hamilton_hill','dickson_lake','disappointment_lake',
           'duffey_lake','gnawed_mountain','highland_valley','mcgillivray_pass',
           'sumallo_river_west','great_bear','upper_squamish','spuzzum_creek','chilliwack_river','tenquille_lake',
           'wahleach_lake','blackwall_peak_pillow')


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

tasmax.file <- paste0('/storage/data/projects/rci/data/winter_sports/course_and_pillow_prism_tasmax_lapse_rates_',type,'.csv')
tasmin.file <- paste0('/storage/data/projects/rci/data/winter_sports/course_and_pillow_prism_tasmin_lapse_rates_',type,'.csv')
pr.file <- paste0('/storage/data/projects/rci/data/winter_sports/course_and_pillow_prism_pr_lapse_rates_',type,'.csv')

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
