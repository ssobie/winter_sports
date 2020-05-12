##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
library(scales)
##---------------------------------------------------------------------------------
##Mapping Component

make_van_whistler_plot <- function(var.name,plot.type,var.title,plot.data,plot.file,shared.range=NULL,mark) {
   if (!is.null(shared.range)) {
      map.range <- shared.range
   } else {
      map.range <- range(as.matrix(plot.data),na.rm=T)
   }

   class.breaks <- c(0,25,50,100,150,200,250,300,400,500,1000)
   ##get.class.breaks(var.name,plot.type,map.range,manual.breaks='')
   colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=map.range,
                                       my.bp=0,class.breaks=class.breaks,
                                       type=plot.type)
   map.class.breaks.labels <- get.class.break.labels(class.breaks)
   white.class.breaks <- c(class.breaks,100000)                                       
   white.colour.ramp <- c(colour.ramp,'#FFFFFF')


   alb.crs <- "+init=epsg:4326"
   lons <- c(-124.0,-123.5, -123.0,-122.5,-122.0,-121.5,-121.0,-120.5,-120.0)
   lats <- c(  48.5,  49.0,  49.5,  50.0,  50.5,  51.0,  51.5,  52.0)

   shape.dir <- '/storage/data/projects/rci/data/assessments/shapefiles/bc'
   bc.shp <- get.region.shape('bc',shape.dir)

   map.extent <- extent(c(-123.5,-121.330,48.9,50.396))
   crop.extent <- map.extent

   plot.window.xlim <- c(map.extent@xmin,map.extent@xmax)
   plot.window.ylim <- c(map.extent@ymin,map.extent@ymax)

   ##par(mar=c(6,6,5,1))
   ##par(mar=mark)
   plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
   bg='white',# 'gray94',
   xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main='',##plot.title,
   cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mgp=c(3.5,2,0),axes=F)
   rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightgray')

   axis(2,at=lats,label=lats,cex.axis=2.25)
   axis(1,at=lons,label=lons,cex.axis=2.25)

   bc.overlay <- 'north_america_state_provincial_boundaries'
   borders.shp <- readOGR(shape.dir, bc.overlay, stringsAsFactors=F, verbose=F)
   wco.shp <- readOGR('/storage/data/projects/rci/data/assessments/shapefiles/bc_common',
                           'ocean_mask', stringsAsFactors=F, verbose=F)
   shape.dir <- paste0('/storage/data/projects/rci/data/assessments/metro_van/shapefiles/')
   region <- 'metro_van'
   region.shp <- spTransform(get.region.shape(region,shape.dir),CRS("+init=epsg:4326")) 

   shade.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'
   shade <- raster(paste0(shade.dir,'van_whistler_hillshade.tif'))


   glacier.dir <- '/storage/data/gis/basedata/randolph_glacier_inventory/v32/02_rgi32_WesternCanadaUS'
   glacier.name <- '02_rgi32_WesternCanadaUS'
   glacier.shp <- spTransform(get.region.shape(glacier.name,glacier.dir),CRS("+init=epsg:4326"))

   image(shade,add=T,col = grey(1:100/100))

   image(plot.data,add=T,col=alpha(white.colour.ramp,0.8),breaks=white.class.breaks)

   plot(spTransform(borders.shp,CRS(alb.crs)),add=TRUE,border='black',cex=0.5)
   plot(spTransform(glacier.shp,CRS(alb.crs)),add=TRUE,col=alpha('white',0.85),border=alpha('white',0.85),cex=0.5)
   plot(spTransform(region.shp,CRS(alb.crs)),add=TRUE,border='black',cex=0.5)
   plot(spTransform(wco.shp,CRS(alb.crs)),add=TRUE,col='gray',lwd=0.5)
      
##   abline(v=lons,lty=3,col='gray',lwd=0.7)
##   abline(h=lats,lty=3,col='gray',lwd=0.7)

   ##par(xpd=NA)
   legend('topright', col = "black", legend=rev(map.class.breaks.labels), pch=22, pt.bg = rev(alpha(colour.ramp,0.8)),
         pt.cex=1.55, y.intersp=0.8, title.adj=0.2, title='mm', xjust=0, cex=1.8)

   box(which='plot',lwd=2)

##   dev.off()
##browser()
}




##-------------------------------------------------------------------------

get.region.shape <- function(region,shape.dir) {
  region.shp <- readOGR(shape.dir, region, stringsAsFactors=F, verbose=F)
  return(region.shp)
}


extract_climatology_data <- function(var.name,interval,clim,
                                     read.dir,tmp.dir,gcm.list,ix) {

   gcm.ens <- c()

   ##Create Ensemble
   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      print(gcm)
      gcm.dir <- paste0(read.dir,gcm,'/')
      var.files <- list.files(path=gcm.dir,pattern=paste0(var.name,'_',clim,'_climatology'))
      gcm.file <- var.files[grep(interval,var.files)]
      print('GCM')
      print(gcm.file)
      file.copy(from=paste0(gcm.dir,gcm.file),to=tmp.dir,overwrite=TRUE)
      gcm.raster <- subset(brick(paste0(tmp.dir,gcm.file)),ix)

      if (is.null(gcm.ens)) {
         gcm.ens <- gcm.raster
      } else {
         gcm.ens <- stack(gcm.ens,gcm.raster)
      }
   }
   rv <- gcm.ens
   return(rv)
}


##-------------------------------------------------------------------------
##Snow model GCMs
tmp.dir <- '/local_temp/ssobie/snow/'
dir.create(tmp.dir,recursive=T)
gcm.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/'
gcm.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')
time <- 'annual'
ix <- 1
past.ens <- extract_climatology_data(var.name='swe',interval='1971-2000',
                                       clim=paste0(time,'_maximum'),
                                       read.dir=gcm.dir,tmp.dir=tmp.dir,
                                       gcm.list=gcm.list,ix=ix)
past.ens.avg <- calc(past.ens,mean)

##-------------------------------------------------------------------------


##---------------------------------------------------------------------------------

##Comparison Plot
if (1==1)  {
plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/ensemble.snow.test.png')
plot.title <- 'ENS Snow'

  png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
  ##par(mfrow=c(2,3))
  ##par(mar=c(0,0,0,0),oma=c(6,6,4,4))
  ##par(mgp=c(4,1.5,0))
  plot.title <- 'ERA5 Total Hazard'

  make_van_whistler_plot(var.name='snowdepth',plot.type='past',plot.title,past.ens.avg,past.plot.file,
                      mark=c(0,6,5,0))
  dev.off()
}

browser()

##---------------------------------------------------------------------------------
##Valid Observations Plot
if (1==0) {
plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/modis.valid.data.png'
plot.title <- 'MODIS Observable Days'
map.range <- range(cover.valid,na.rm=T)
leg.title <- 'Days'
class.breaks <- c(0,1000,seq(1250,2500,by=250)) ##get.class.breaks('tasmax',type='past',map.range,manual.breaks='')
map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=FALSE)
colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)
vw.plot(valid.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=FALSE)
}

##---------------------------------------------------------------------------------
##Percent Difference Sum
if (1==0) {
plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/modis-',model,'.abs.diff.png')
plot.title <- paste0(model,'-MODIS Percent Diff Sum')
leg.title <- 'Diff. (%)'
snow.prct <- (model.snow-modis.snow)/modis.snow*100
sp.raster <- list(x=lon,y=lat,z=snow.prct)
map.range <- range(snow.prct,na.rm=T)
class.breaks <- c(-5000,-100,-80,-60,-40,-20,-10,0,10,20,50,100,5000) ##c(-1000,-500,-250,-100,-50,-25,0,25,50,100,200,300,1000)  ##
map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=TRUE,lesser.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)
vw.plot(sp.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=FALSE)
}

##---------------------------------------------------------------------------------
##Absolute Difference Sum
if (1==0) {
plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/modis-era.absolute.diff.png'
plot.title <- paste0(model,'-MODIS Absolute Diff Sum')
leg.title <- 'Diff. (Days)'
snow.prct <- (model.snow-modis.snow)##/modis.snow*100
sp.raster <- list(x=lon,y=lat,z=snow.prct)
map.range <- range(snow.prct,na.rm=T)
class.breaks <- c(-5000,-100,-50,-20,0,20,50,100,200,300,5000) ##get.class.breaks('pr',type='anom',map.range,manual.breaks='')
map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=TRUE,lesser.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)
vw.plot(sp.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE)
}

##---------------------------------------------------------------------------------
##MODIS Total Snow Cover
if (1==0) {
plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/modis.total.snow.cover.days.png'
plot.title <- 'MODIS Total Snow Cover Days'
leg.title <- 'Days'
map.range <- range(modis.snow,na.rm=T)
class.breaks <- c(0,10,20,50,100,300,900,1500,5000) ##get.class.breaks('pr',type='prct',map.range,manual.breaks='')
map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)
vw.plot(modis.snow.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE)
}

##---------------------------------------------------------------------------------
##Four panel figure
if (1==0) {
plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/modis.metro.van.comparison.data.2018.png'
png(plot.file,width=1800,height=1800)
par(mfrow=c(3,2))

##Observable days
plot.title <- 'MODIS Observable Days' 
load(paste0(save.dir,'modis.observable.days.2018.RData'))
map.range <- c(0,3500) ##range(cover.valid,na.rm=T)
leg.title <- 'Days'
class.breaks <- c(0,1000,seq(1250,3500,by=250)) ##get.class.breaks('tasmax',type='past',map.range,manual.breaks='')
map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=FALSE)
colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)
vw.plot(valid.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=FALSE)

##Total Snow Cover
plot.title <- 'MODIS Days with Snow Cover'
leg.title <- 'Days'
load(paste0(save.dir,'modis.snow.days.2018.RData'))
map.range <- c(0,1500) ##range(modis.snow,na.rm=T)
class.breaks <- c(0,10,20,50,100,300,900,1500,5000) ##get.class.breaks('pr',type='prct',map.range,manual.breaks='')
map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)
vw.plot(modis.snow.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=FALSE)

##Success Rate for ERA
plot.title <- paste0('ERA-Interim - MODIS Snow Cover Comparison')
load(paste0(save.dir,'ERA.success.rate.2018.RData'))
map.range <- c(70,100)
leg.title <- 'Success (%)'
class.breaks <- c(0,get.class.breaks('snc',type='past',map.range,manual.breaks=''))
map.class.breaks.labels <- get.class.break.labels(class.breaks,lesser.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='snd',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)

vw.plot(ratio.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE)

##Success Rate for NCEP2
plot.title <- paste0('NCEP2 - MODIS Snow Cover Comparison')
load(paste0(save.dir,'NCEP2.success.rate.2018.RData'))
map.range <- c(70,100)
leg.title <- 'Success (%)'
class.breaks <- c(0,get.class.breaks('snc',type='past',map.range,manual.breaks=''))
map.class.breaks.labels <- get.class.break.labels(class.breaks,lesser.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='snd',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)

vw.plot(ratio.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE)


plot.title <- paste0('ERA-Interim - MODIS Percent Diff Sum')
load(paste0(save.dir,'ERA.snow.days.2018.RData'))
model.snow <- model.snow.raster$z
modis.snow <- modis.snow.raster$z
lon <- model.snow.raster$x
lat <- model.snow.raster$y
leg.title <- 'Diff. (%)'
snow.prct <- (model.snow-modis.snow)/modis.snow*100
sp.raster <- list(x=lon,y=lat,z=snow.prct)
map.range <- range(snow.prct,na.rm=T)
class.breaks <- c(-5000,-100,-80,-60,-40,-20,-10,0,10,20,50,100,5000) ##c(-1000,-500,-250,-100,-50,-25,0,25,50,100,200,300,1000)  ##
map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=TRUE,lesser.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)
vw.plot(sp.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE)
}


##plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/modis-NCEP2.prct.diff.png'
##png(plot.file,width=1600,height=1600)

plot.title <- paste0('NCEP2 - MODIS Percent Diff Sum')
load(paste0(save.dir,'NCEP2.snow.days.2018.RData'))
model.snow <- model.snow.raster$z
modis.snow <- modis.snow.raster$z
lon <- model.snow.raster$x
lat <- model.snow.raster$y

leg.title <- 'Diff. (%)'
snow.prct <- (model.snow-modis.snow)/modis.snow*100
sp.raster <- list(x=lon,y=lat,z=snow.prct)
map.range <- range(snow.prct,na.rm=T)
class.breaks <- c(-5000,-100,-80,-60,-40,-20,-10,0,10,20,50,100,5000) ##c(-1000,-500,-250,-100,-50,-25,0,25,50,100,200,300,1000)  ##
map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=TRUE,lesser.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)
vw.plot(sp.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE)




dev.off()