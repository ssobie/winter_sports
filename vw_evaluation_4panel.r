##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/ncc_panel_map_support.r',chdir=T)

library(raster)

##-------------------------------------------------------------------------------

extract_climatology_data <- function(var.name,interval,clim,
                                     read.dir,clim.dir,gcm.list) {

   gcm.ens <- c()

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

      gcm.raster <- subset(brick(paste0(gcm.dir,gcm.file)),1)

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
base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/'

gcm.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')

var.name <- 'swe'
clim <- 'seasonal_mean'

snow.file <- paste0(base.dir,'calibrated_PNWNAmet_PNWNAmet_prism_tps/',
                    'swe_seasonal_mean_climatology_BCCAQ2-PRISM_PNWNAmet_PNWNAmet_1950-2012.nc')
snow.seas <- brick(snow.file)*1000
pnw.1951.2012 <- subset(snow.seas,1)

snow.file <- paste0(base.dir,'calibrated_ERA5_PNWNAmet_prism_tps/',
                    'swe_seasonal_mean_climatology_BCCAQ2-PRISM_ERA5_PNWNAmet_1981-2018.nc')
snow.seas <- brick(snow.file)*1000
era.1981.2018 <- subset(snow.seas,1)

past.ens <- extract_climatology_data(var.name=var.name,interval='1950-2012',
                                     clim=clim,
                                     read.dir=base.dir,clim.dir='standard_climatologies',
                                     gcm.list=gcm.list)

##past.ens.avg <- calc(past.ens,mean)*1000
past.ens.avg <- calc(past.ens,fun=function(x){quantile(x,0.1,na.rm=T,names=FALSE)})*1000

percent.diff <- (past.ens.avg-pnw.1951.2012)/pnw.1951.2012*100
percent.diff[is.infinite(percent.diff)] <- NA
percent.diff[percent.diff > 10 ] <- 10

snow.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/')
site.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'
active.courses <- read.csv(paste0(site.dir,'active_snow_courses.csv'),header=T,as.is=T)
inactive.courses <- read.csv(paste0(site.dir,'inactive_snow_courses.csv'),header=T,as.is=T)
snow.pillows <- read.csv(paste0(site.dir,'snow_pillow_locations.csv'),header=T,as.is=T)

##---------------------------------------------------

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/swe.winter.mean.evaluation.4.panel.10th.percentile.png')
##plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/swe.spring.mean.evaluation.4.panel.png')
plot.title <- '' ##'Annual Peak SWE'
##class.breaks <- c(0,25,50,100,150,200,250,300,400,500,1000,1500,2000)
class.breaks <- c(0,25,50,100,150,200,250,500,1000,1500,2000,2500,3000,3500)
percent.breaks <- seq(-100,10,10)

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
par(mfrow=c(2,2))
par(mar=c(0,0,0,0),oma=c(6,6,4,14))

rv.swe <- make_van_whistler_panel_plot(var.name='swe',plot.type='past',plot.title,pnw.1951.2012,
                       class.breaks=class.breaks,y.axis=TRUE,letter='A')
rv <- make_van_whistler_panel_plot(var.name='swe',plot.type='past',plot.title,
                       era.1981.2018,
                       class.breaks=class.breaks,letter='B')
par(xpd=NA)
legend('topright',inset=c(-0.29,0), col = "black",
        legend=rv.swe$labels, pch=22, pt.bg = rv.swe$cols,
        pt.cex=3.55, y.intersp=0.8, title.adj=0.2, title='SWE (mm)', xjust=0, cex=1.8,box.lwd=2)
par(xpd=FALSE)

rv <- make_van_whistler_panel_plot(var.name='swe',plot.type='past',plot.title,
                       past.ens.avg,
                       class.breaks=class.breaks,x.axis=TRUE,y.axis=TRUE,letter='C')
rv <- make_van_whistler_panel_plot(var.name='swe',plot.type='diff',plot.title,
                       percent.diff,add.legend=TRUE,leg.title='%',
                       class.breaks=percent.breaks,x.axis=TRUE,letter='D')

mtext("Longitude (\u00B0E)",side=1,outer=TRUE,cex=2.25,line=4.6)
mtext("Latitude (\u00B0N)",side=2,outer=TRUE,cex=2.25,line=3.6)

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

dev.off()
