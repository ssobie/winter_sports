##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/ncc_map_support.r',chdir=T)
 

make_cal_map <- function(var.name,plot.data,leg.title,class.breaks=NULL,mark) {

    map.range <- range(class.breaks)
    colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                        my.bp=0,class.breaks=class.breaks,
                                        type='past')
   map.class.breaks.labels <- get.class.break.labels(class.breaks)

   alb.crs <- "+init=epsg:4326"
   lons <- c(-124.0,-123.5, -123.0,-122.5,-122.0,-121.5,-121.0,-120.5,-120.0)
   lats <- c(48.75, 49.0, 49.25, 49.5, 49.75, 50.0, 50.25, 50.5, 50.75)

   map.extent <- extent(c(-123.7,-120.65,48.9,50.8))
   plot.window.xlim <- c(map.extent@xmin,map.extent@xmax)
   plot.window.ylim <- c(map.extent@ymin,map.extent@ymax)

   plot.title <- ''
   if (mark[2]==6 & mark[3]==5){plot.title=='PNWNAmet'}
   if (mark[3]==5 & mark[4]==1){plot.title=='ERA5'}
   
   plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
   bg='white',# 'gray94',
   xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main=plot.title,
   cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mgp=c(3.5,2,0),axes=F)
   rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightgray')

   if (mark[2] == 6) {
     axis(2,at=lats,label=lats,cex.axis=3.25,mgp=c(4,2.5,0))
   }

   if (mark[1]==5) {
     axis(1,at=lons,label=lons,cex.axis=3.25,mgp=c(4,3,0))
   }

   shape.dir <- '/storage/data/projects/rci/data/assessments/shapefiles/bc'
   bc.overlay <- 'north_america_state_provincial_boundaries'
   borders.shp <- readOGR(shape.dir, bc.overlay, stringsAsFactors=F, verbose=F)
   wco.shp <- readOGR('/storage/data/projects/rci/data/assessments/shapefiles/bc_common',
                           'ocean_mask', stringsAsFactors=F, verbose=F)
   shape.dir <- paste0('/storage/data/projects/rci/data/assessments/metro_van/shapefiles/')
   region <- 'metro_van'
   region.shp <- spTransform(get.region.shape(region,shape.dir),CRS("+init=epsg:4326"))
   whistler.dir <- '/storage/data/projects/rci/data/assessments/whistler/shapefiles/'
   whistler.shp  <- spTransform(get.region.shape('WhistlerLandscapeUnit',whistler.dir),CRS("+init=epsg:4326"))
   shade.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'
   shade <- raster(paste0(shade.dir,'van_whistler_hillshade.tif'))

   wash.shp <- spTransform(readOGR(shade.dir, 'washington', stringsAsFactors=F, verbose=F),CRS("+init=epsg:4326"))

   image(shade,add=T,col = grey(1:100/100))
   image(plot.data,add=T,col=alpha(colour.ramp,0.8),breaks=class.breaks)
   plot(borders.shp,add=T)
   site.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'
   active.courses <- read.csv(paste0(site.dir,'active_snow_courses.csv'),header=T,as.is=T)
   inactive.courses <- read.csv(paste0(site.dir,'inactive_snow_courses.csv'),header=T,as.is=T)
   snow.pillows <- read.csv(paste0(site.dir,'snow_pillow_locations.csv'),header=T,as.is=T)

   points(active.courses$Lon,active.courses$Lat,pch=24,col='black',bg='black',cex=2.5)
   points(inactive.courses$Lon,inactive.courses$Lat,pch=25,col='black',bg='black',cex=2.5)
   points(snow.pillows$Lon,snow.pillows$Lat,pch=23,col='black',bg='black',cex=2.5)
   
   if (leg.title=='Scale' & mark[4]==1) {inset <- c(-0.22,0)}
   if (leg.title=='Slope' & mark[4]==1) {inset <- c(-0.24,0)}
   if (leg.title=='Freq.' & mark[4]==1) {inset <- c(-0.24,0)}

   if (mark[4]==1) {
      par(xpd=NA)
      legend('topright',inset=inset, col = "black", legend=rev(map.class.breaks.labels), pch=22, pt.bg = rev(alpha(colour.ramp,0.8)),
            pt.cex=4, y.intersp=0.8, title.adj=0.2, title=leg.title, xjust=0, cex=2.8,box.lwd=2)
      par(xpd=FALSE)
   }
   box(which='plot',lwd=2)
}

##--------------------------------------------------------------------------------------------

##ERA5 - PRISM
cal.dir <- '/storage/data/projects/rci/data/winter_sports/'
scale.file <- paste0(cal.dir,'scale_hyper_snow_calibrated_parameter_ERA5_prism_TPS.nc')
scale.nc <- nc_open(scale.file)
lon <- ncvar_get(scale.nc,'lon')
lat <- ncvar_get(scale.nc,'lat')
scale.era5 <- ncvar_get(scale.nc,'scale')
nc_close(scale.nc)

scale.era5.raster <-  list(x=lon,y=lat,z=scale.era5)

era5.scale <- brick(scale.file)
era5.slope <- brick(paste0(cal.dir,'slope_hyper_snow_calibrated_parameter_ERA5_prism_TPS.nc'))
era5.freq <- brick(paste0(cal.dir,'freq_hyper_snow_calibrated_parameter_ERA5_prism_TPS.nc'))

pnw.scale <- brick(paste0(cal.dir,'scale_hyper_snow_calibrated_parameter_PNWNAmet_prism_TPS.nc'))
pnw.slope <- brick(paste0(cal.dir,'slope_hyper_snow_calibrated_parameter_PNWNAmet_prism_TPS.nc'))
pnw.freq <- brick(paste0(cal.dir,'freq_hyper_snow_calibrated_parameter_PNWNAmet_prism_TPS.nc'))

plot.file <- paste('/storage/data/projects/rci/data/winter_sports/plots/snow.model.calibration.parameters.TPS.png') 

png(file=plot.file,width=10.5,height=10,units='in',res=600,pointsize=6,bg='white')
par(mfrow=c(3,2))
par(mar=c(0,0,0,0),oma=c(10,9,5,17))

class.breaks <- seq(44,55,1)
make_cal_map(var.name='scale',pnw.scale,'Scale',
                       class.breaks=class.breaks,mark=c(0,6,5,0))
make_cal_map(var.name='scale',era5.scale,'Scale',
                       class.breaks=class.breaks,mark=c(0,0,5,1))
class.breaks <- seq(0.0,2,0.2)
make_cal_map(var.name='slope',pnw.slope,'Slope',
                       class.breaks=class.breaks,mark=c(0,6,0,0))
make_cal_map(var.name='slope',era5.slope,'Slope',
                       class.breaks=class.breaks,mark=c(0,0,0,1))
class.breaks <- seq(-1.5,5.5,0.5)
make_cal_map(var.name='freq',pnw.freq,'Freq.',
                       class.breaks=class.breaks,mark=c(5,6,0,0))
make_cal_map(var.name='freq',era5.freq,'Freq.',
                       class.breaks=class.breaks,mark=c(5,0,0,1))


mtext("Longitude (\u00B0E)",side=1,outer=TRUE,cex=2.6,line=6.6)
mtext("Latitude (\u00B0N)",side=2,outer=TRUE,cex=2.6,line=5.6)


dev.off()

browser()


