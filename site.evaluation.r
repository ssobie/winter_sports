##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)
library(abind)

source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##-------------------------------------------------------------------------------

extract_climatology_data <- function(var.name,interval,clim,
                                     read.dir,clim.dir,gcm.list) {

   gcm.ens <- array(NA,c(480,323,length(gcm.list)))

   ##Create Ensemble
   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      print(gcm)
      gcm.dir <- paste0(read.dir,'calibrated_',gcm,'_PNWNAmet_prism_tps/',clim.dir,'/')
      var.files <- list.files(path=gcm.dir,pattern=paste0(var.name,'_',clim,'_climatology'))
      if (interval == 'one') {
         gcm.file <- var.files[1]
      } else if (interval == 'two') {
         gcm.file <- var.files[2]
      } else if (interval == 'three') {
         gcm.file <- var.files[3]
      } else {
         gcm.file <- var.files[grep(interval,var.files)]
      }
      print('GCM')
      print(gcm.file)
      nc <- nc_open(paste0(gcm.dir,gcm.file))
      gcm.ens[,,g] <- ncvar_get(nc,var.name,start=c(1,1,1),count=c(-1,-1,1))

      nc_close(nc)
   }
   rv <- gcm.ens*1000
   return(rv)
}
##--------------------------------------------------------------

get_diff_values <- function(sites,asps,gcm.list,lon,lat,pnw.snow,era5.snow,gcm.pnw,gcm.era5) {

   course.site.pnw <- course.site.era5 <- matrix(NA,nrow=length(sites),ncol=length(gcm.list))
   pillow.site.pnw <- pillow.site.era5 <- matrix(NA,nrow=length(sites),ncol=length(gcm.list))

   ##Coordinates to match site coords

   for (i in seq_along(sites)) {

       site <- sites[i]
       print(site)
       coords <- get_coordinates(site)
       lon.ix <- which.min(abs(coords[1]-lon))
       lat.ix <- which.min(abs(coords[2]-lat))
       course.site.pnw[i,] <- (as.numeric(gcm.pnw[lon.ix,lat.ix,]) - as.numeric(pnw.snow[lon.ix,lat.ix]))
       course.site.era5[i,] <- (as.numeric(gcm.era5[lon.ix,lat.ix,]) - as.numeric(era5.snow[lon.ix,lat.ix]))

   }        

   for (i in seq_along(asps)) {
      asp <- asps[i]
      coords <- get_coordinates(asp)
      lon.ix <- which.min(abs(coords[1]-lon))
      lat.ix <- which.min(abs(coords[2]-lat))
      pillow.site.pnw[i,] <- (as.numeric(gcm.pnw[lon.ix,lat.ix,]) - as.numeric(pnw.snow[lon.ix,lat.ix]))
      pillow.site.era5[i,] <- (as.numeric(gcm.era5[lon.ix,lat.ix,]) - as.numeric(era5.snow[lon.ix,lat.ix]))
   }        

   all.pnw.sites <- rbind(course.site.pnw,pillow.site.pnw)
   all.era5.sites <- rbind(course.site.era5,pillow.site.era5)

   rv <- list(era=all.era5.sites,pnw=all.pnw.sites)           
   return(rv)
}

##---------------------------------------------------------------------------

##Snow model GCMs
base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/'

gcm.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')

var.name <- 'swe'
clim <- 'seasonal_mean'

snow.file <- paste0(base.dir,'calibrated_PNWNAmet_PNWNAmet_prism_tps/',
                    'swe_seasonal_mean_climatology_BCCAQ2-PRISM_PNWNAmet_PNWNAmet_1950-2012.nc')

nc <- nc_open(snow.file)
lon <- ncvar_get(nc,'lon')
lat <- ncvar_get(nc,'lat')
pnw.1951.2012 <- ncvar_get(nc,var.name,start=c(1,1,1),count=c(-1,-1,1))*1000
nc_close(nc)


snow.file <- paste0(base.dir,'calibrated_ERA5_PNWNAmet_prism_tps/',
                    'swe_seasonal_mean_climatology_BCCAQ2-PRISM_ERA5_PNWNAmet_1981-2018.nc')
nc <- nc_open(snow.file)
era5.1981.2018 <- ncvar_get(nc,var.name,start=c(1,1,1),count=c(-1,-1,1))*1000
nc_close(nc)

past.ens.1950 <- extract_climatology_data(var.name=var.name,interval='1950-2012',
                                     clim=clim,
                                     read.dir=base.dir,clim.dir='standard_climatologies',
                                     gcm.list=gcm.list)

past.ens.1980 <- extract_climatology_data(var.name=var.name,interval='1981-2018',
                                     clim=clim,
                                     read.dir=base.dir,clim.dir='standard_climatologies',
                                     gcm.list=gcm.list)
 

##*******************************************************************************


sites <- c('shovelnose_mountain','brookmere','lightning_lake','callaghan','orchid_lake',
           'palisade_lake','grouse_mountain','dog_mountain','stave_lake','nahatlatch',
           'wahleach','klesilkwa','hamilton_hill','dickson_lake','disappointment_lake',
           'duffey_lake','gnawed_mountain','highland_valley','mcgillivray_pass',
           'sumallo_river_west','great_bear')
ilen <- length(sites)
asps <- c('upper_squamish','spuzzum_creek','chilliwack_river','tenquille_lake',
           'wahleach_lake','blackwall_peak_pillow')

site.names <- c('Shovelnose', ##\nMountain',
                'Brookmere',
                'Lightning', ##\nLake',
                'Callaghan',
                'Orchid', ##\nLake',
                'Palisade', ##\nLake',
                'Grouse', ##\nMountain',
                'Dog', ##\nMountain',
                'Stave', ##\nLake',
                'Nahatlatch',
                'Wahleach',
                'Klesilkwa',
                'Hamilton', ##\nHill',
                'Dickson', ##\nLake',
                'Disappoint.', ##ment', ##\nLake',
                'Duffey', ##\nLake',
                'Gnawed', ##\nMountain',
                'Highland', ##\nValley',
                'Mcgillivray', ##\nPass',
                'Sumallo', ##\nWest',          
                'Bear', ##(Great)
                'Squamish', ##Upper\n
                'Spuzzum', ##\nCreek',
                'Chilliwack', ##\nRiver',
                'Tenquille', ##\nLake',
                'Wahleach L.', ##\nLake',
                'Blackwall') ##\nPeak')

course.lons <- rep(0,length(sites))
pillow.lons <- rep(0,length(asps))
course.lats <- rep(0,length(sites))
pillow.lats <- rep(0,length(asps))

course.elevs <- rep(0,length(sites))
pillow.elevs <- rep(0,length(asps))

for (i in seq_along(sites)) {
    site <- sites[i]
    coords <- get_coordinates(site)
    course.lons[i] <- coords[1]
    course.lats[i] <- coords[2]
    course.elevs[i] <- coords[3]
}

for (i in seq_along(asps)) {
    asp <- asps[i]
    coords <- get_coordinates(asp)
    pillow.lons[i] <- coords[1]
    pillow.lats[i] <- coords[2]
    pillow.elevs[i] <- coords[3]
}

lons <- c(course.lons,pillow.lons)
lats <- c(course.lats,pillow.lats)
elevs <- c(course.elevs,pillow.elevs)
alen <- length(sites)+length(asps)

cal.diff <- get_diff_values(sites,asps,gcm.list,lon,lat,pnw.1951.2012,era5.1981.2018,
                            past.ens.1950,past.ens.1980)


##-------------------------------------------------------------------

##Loop over sites
plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'

leg.title <- 'mm'
ranked.lons <- order(lons)
ranked.elevs <- order(elevs)

png(file=paste0(plot.dir,'pnwnamet.era5.gcm.evaluation.2020.png'),width=10,height=6,units='in',res=600,pointsize=6,bg='white')

layout(mat = matrix(c(1,2),
                    nrow = 2,
                    ncol = 1),
             heights = c(1.6, 1.4),    # Heights of the two rows
             widths = 1)     # Widths of the two columns

par(mar=c(9,5,2,3))
plot(0:alen,0:alen,xlab='',ylab='SWE Bias (mm)',yaxs='i',
     col='white',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,alen),ylim=c(-600,400),axes=FALSE)
axis(1,at=1:alen,site.names[ranked.lons],cex=1.75,cex.axis=1.75,las=2)
axis(2,at=c(-400,-200,200,400),label=c(-400,-200,200,400),cex=1.75,cex.axis=1.75)
abline(h=seq(-400,400,200),lty=2,col='gray',lwd=2)
abline(v=1:alen,col='gray')
abline(h=0,lwd=1.5)
for (j in 1:alen) {
    print(lons[ranked.lons[j]])
    boxplot(at=j-0.175,x=cal.diff$era[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='blue',border='black')
    boxplot(at=j+0.175,x=cal.diff$pnw[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='green',border='black')
}
##text(x=16.6,y=550,'SNODAS',cex=1.25)
##legend('bottomright',leg=c('ERA5','PNWNAmet','SNODAS'),col=c('blue','green','red'),pch=15,cex=1.5)
box(which='plot',lwd=1.5)

sites <- c('blackwall_peak_course','boston_bar_lower','boston_bar_upper','burwell_lake',
           'chapman_creek','cornwall_hills','diamond_head','edwards_lake',
           'hollyburn','hope',
           'loch_lomond','lytton','mount_seymour','new_tashme',
           'ottomite','pavilion_mountain',
           'sumallo_river','tenquille_course','whistler_mountain','wolverine_creek' )

elen <- length(sites)

site.names <- c('Blackwall',
                'Boston U.',
                'Boston L.',
                'Burwell',
                'Chapman', 
                'Cornwall',
                'Diamond', 
                'Edwards', 
                'Hollyburn',
                'Hope',
                'Lomond',
                'Lytton', 
                'Seymour',
                'Tashme',
                'Ottomite',
                'Pavilion',
                'Sumallo',
                'Tenquille',
                'Whistler', 
                'Wolverine')

course.lons <- rep(0,length(sites))
course.lats <- rep(0,length(sites))

course.elevs <- rep(0,length(sites))

for (i in seq_along(sites)) {
    site <- sites[i]
    coords <- get_coordinates(site)
    course.lons[i] <- coords[1]
    course.lats[i] <- coords[2]
    course.elevs[i] <- coords[3]
}

lons <- course.lons
lats <- course.lats
elevs <- course.elevs
elen <- length(sites)
ranked.lons <- order(lons)
ranked.elevs <- order(elevs)

val.diff <- get_diff_values(sites,asps,gcm.list,lon,lat,pnw.1951.2012,era5.1981.2018,
                            past.ens.1950,past.ens.1980)

##slen <- length(sites)
##val.mae <- lapply(val.mae,function(x,slen){x[1:slen]},slen)

par(mar=c(9,5,2,9))
plot(0:elen,0:elen,xlab='',ylab='SWE Bias (mm)',yaxs='i',
     col='white',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,elen),ylim=c(-600,600),axes=FALSE)
axis(1,at=1:elen,site.names[ranked.lons],cex=1.75,cex.axis=1.75,las=2)
axis(2,at=c(-400,-200,200,400),label=c(-400,-200,200,400),cex=1.75,cex.axis=1.75)
abline(h=seq(-400,400,200),lty=2,col='gray',lwd=2)
abline(v=1:alen,col='gray')
abline(h=0,lwd=1.5)
for (j in 1:elen) {
    print(lons[ranked.lons[j]])
    boxplot(at=j-0.175,x=val.diff$era[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='blue',border='black') 
    boxplot(at=j+0.175,x=val.diff$pnw[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='green',border='black') 
}
par(xpd=NA)
legend('topright',inset=c(-0.08,0),leg=c('ERA5','PNW'),col=c('blue','green'),pch=c(15,15),pt.cex=c(2,2),cex=1.5,box.lwd=1.5)
par(xpd=FALSE)

box(which='plot',lwd=1.5)


dev.off()

browser()

png(filename=paste0(plot.dir,'era.ncep2.all.sites.swe.mae.2018.png'),width=1000,height=700)
layout(mat = matrix(c(1,2), 
                    nrow = 2, 
                    ncol = 1),
             heights = c(1, 1.5),    # Heights of the two rows
             widths = 1)     # Widths of the two columns

##par(mfrow=c(2,1))

par(mar=c(0,5,3,3))
plot(0:alen,0:alen,xlab='',ylab='MAE (mm)',yaxs='i',
     col='white',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,alen),ylim=c(0,600),axes=FALSE)
##axis(1,at=1:alen,site.names[ranked.lons],cex=1.75,cex.axis=1.75,las=2)
axis(2,at=c(0,200,400,600),c(0,200,400,600),cex=1.75,cex.axis=1.75)
abline(h=seq(0,1000,100),lty=2,col='gray',lwd=2)
abline(v=1:alen,col='gray')
for (j in 1:alen) {
    print(lons[ranked.lons[j]])
    boxplot(at=j-0.175,x=era.snodas.mae[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='blue',border='blue')
    boxplot(at=j+0.175,x=ncep2.snodas.mae[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='green',border='green')
}

text(x=16.6,y=550,'SNODAS',cex=1.25)
legend('topleft',leg=c('ERA-I','NCEP2'),col=c('blue','green'),pch=15,cex=1.5)

box(which='plot')

par(mar=c(10,5,0.1,3))
plot(0:alen,0:alen,xlab='',ylab='MAE (mm)',yaxs='i',
     col='white',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,alen),ylim=c(0,1200),axes=FALSE)
axis(1,at=1:alen,site.names[ranked.lons],cex=1.75,cex.axis=1.75,las=2)
axis(2,at=seq(0,1200,200),seq(0,1200,200),cex=1.75,cex.axis=1.75)
abline(h=seq(0,2000,200),lty=2,col='gray',lwd=2)
abline(v=1:alen,col='gray')
for (j in 1:alen) {
    print(lons[ranked.lons[j]])
    boxplot(at=j-0.175,x=all.era.mae[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='blue',border='blue')
    boxplot(at=j+0.175,x=all.ncep2.mae[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='green',border='green')
    points(x=j,y=all.snodas.mae[ranked.lons[j]],pch='-',cex=5,col='red')
}
text(x=16.0,y=1100,'Courses and Pillows',cex=1.25)
legend('topleft',leg=c('ERA-I','SNODAS','NCEP2'),col=c('blue','red','green'),pch=15,cex=1.5)

box(which='plot')

dev.off()
