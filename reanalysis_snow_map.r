##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/ncc_map_support.r',chdir=T)

library(raster)

##-------------------------------------------------------------------------------

extract_climatology_data <- function(var.name,interval,clim,
                                     read.dir,gcm.list,ix) {

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
      gcm.raster <- subset(brick(paste0(gcm.dir,gcm.file)),ix)

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

site.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'
active.courses <- read.csv(paste0(site.dir,'active_snow_courses.csv'),header=T,as.is=T)
inactive.courses <- read.csv(paste0(site.dir,'inactive_snow_courses.csv'),header=T,as.is=T)
snow.pillows <- read.csv(paste0(site.dir,'snow_pillow_locations.csv'),header=T,as.is=T)

##---------------------------------------------------
model <- 'PNWNAmet'
type <- 'annual_maximum'
interval <- '1980-2012'
seas.ix <- 1

snow.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/calibrated_',model,'_prism_tps/')
snow.file <- paste0(snow.dir,'swe_',type,'_climatology_BCCAQ2-PRISM_',model,'_',interval,'.nc')
snow.seas <- brick(snow.file)*1000
snow.data <- subset(snow.seas,seas.ix)

if (1==1) {

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/',
                     'swe.',type,'.',model,'.courses.map.',interval,'.png')
plot.title <- 'Annual Peak SWE'
##class.breaks <- c(0,25,50,100,150,200,250,300,400,500,1000,1500,2000)
class.breaks <- c(0,25,50,100,150,200,250,500,1000,1500,2000,2500,3000)

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='swe',plot.type='past',plot.title,snow.data,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')

points(active.courses$Lon,active.courses$Lat,pch=24,col='black',bg='green',cex=1.5)
points(inactive.courses$Lon,inactive.courses$Lat,pch=25,col='black',bg='red',cex=1.5)
points(snow.pillows$Lon,snow.pillows$Lat,pch=23,col='black',bg='white',cex=1.5)

dev.off()


browser()
##---------------------------------------------------
##Model Past Peak SWE

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/swe.peak.courses.ENS.1971-2000.png')
plot.title <- 'Peak SWE (1971-2000)'

class.breaks <- c(0,25,50,100,150,200,250,500,1000,1500,2000,2500,3000)

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='swe',plot.type='past',plot.title,past.ens.avg,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')

points(course.left$Lon,course.left$Lat,pch=18,cex=3.5)
points(course.right$Lon,course.right$Lat,pch=18,cex=3.5)
points(snow.pillows$Lon,snow.pillows$Lat,pch=17,cex=2.5)

points(course.left$Lon,course.left$Lat,pch=18,cex=2.0,col='red')
points(course.right$Lon,course.right$Lat,pch=18,cex=2.0,col='red')
points(snow.pillows$Lon,snow.pillows$Lat,pch=17,cex=1.0,col='red')


dev.off()


##---------------------------------------------------
##Model Future Peak SWE

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/swe.peak.courses.ENS.2041-2070.png')
plot.title <- 'Peak SWE (2041-2070)'
##class.breaks <- c(0,25,50,75,100,150,200,250,300,350,400,500)
class.breaks <- c(0,25,50,100,150,200,250,500,1000,1500,2000,2500,3000)

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='swe',plot.type='past',plot.title,proj.ens.avg,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')

points(course.left$Lon,course.left$Lat,pch=18,cex=3.5)
points(course.right$Lon,course.right$Lat,pch=18,cex=3.5)
points(snow.pillows$Lon,snow.pillows$Lat,pch=17,cex=2.5)

points(course.left$Lon,course.left$Lat,pch=18,cex=2.0,col='red')
points(course.right$Lon,course.right$Lat,pch=18,cex=2.0,col='red')
points(snow.pillows$Lon,snow.pillows$Lat,pch=17,cex=1.0,col='red')

dev.off()

}

##---------------------------------------------------
##Model Anomaly Peak SWE

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/swe.peak.courses.anomalies.ENS.2041-2070_1971-2000.png')
plot.title <- 'Peak SWE Anomalies (2041-2070 - 1971-2000)'

class.breaks <- seq(-90,0,10) ##class.breaks <- c(-3000,-2500,-2000,-1500,-1000,-500,-250,-200,-150,-100,-50,-25,0)

proj.ens.avg[proj.ens.avg > 3000] <- 3000
past.ens.avg[past.ens.avg > 3000] <- 3000

diff.ens.avg <- proj.ens.avg - past.ens.avg
prct.ens.avg <- (proj.ens.avg - past.ens.avg)/past.ens.avg *100

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='swe',plot.type='percent',plot.title,prct.ens.avg,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='%')

points(course.left$Lon,course.left$Lat,pch=18,cex=3.5)
points(course.right$Lon,course.right$Lat,pch=18,cex=3.5)
points(snow.pillows$Lon,snow.pillows$Lat,pch=17,cex=2.5)

points(course.left$Lon,course.left$Lat,pch=18,cex=2.0,col='red')
points(course.right$Lon,course.right$Lat,pch=18,cex=2.0,col='red')
points(snow.pillows$Lon,snow.pillows$Lat,pch=17,cex=1.0,col='red')

dev.off()
