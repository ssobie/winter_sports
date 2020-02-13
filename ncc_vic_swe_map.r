##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/ncc_map_support.r',chdir=T)

library(raster)


##-------------------------------------------------------------------------------

extract_vic_climatology_data <- function(var.name,interval,clim,
                                         read.dir,gcm.list,ix) {

   gcm.ens <- c()

   ##Create Ensemble
   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      print(gcm)
      var.files <- list.files(path=read.dir,pattern=paste0(var.name,'_',clim,'_vw_climatology_',gcm))
      gcm.file <- var.files[grep(interval,var.files)]
      print('GCM')
      print(gcm.file)
      gcm.raster <- subset(brick(paste0(read.dir,gcm.file)),ix)
      gcm.raster@file@nodatavalue <- -Inf
##browser()
      if (is.null(gcm.ens)) {
         gcm.ens <- gcm.raster
      } else {
         gcm.ens <- stack(gcm.ens,gcm.raster)
      }
   }
   rv <- gcm.ens
   return(rv)
}


##---------------------------------------------------------------------
##Past simulations
if (1==0) {
vic.dir <- paste0('/storage/data/projects/rci/data/winter_sports/')
vic.file <- paste0(vic.dir,'vic_swe_season_clim.nc')
vic.seas <- brick(vic.file)
vic.winter <- subset(vic.seas,1)
vic.spring <- subset(vic.seas,2)

vic.data <- vic.winter*1000

class.breaks <- c(0,25,50,100,150,200,250,300,400,500,1000,1500,2000)

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/vic.base.swe.map.png')
plot.title <- 'VIC SWE'

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')

make_van_whistler_plot(var.name='swe',plot.type='past',plot.title,vic.data,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')
dev.off()
}

##-----------------------------------------------------------------------
##GCM VIC SWE
read.dir <- paste0('/storage/data/projects/rci/data/winter_sports/VIC/')
time <- 'annual'
ix <- 1
var.name <- 'swe'

gcm.list <- c('CCSM3','CGCM3','CSIRO35','ECHAM5','GFDL2.1','HadCM','HadGEM1','MIROC3.2')

past.ens <- extract_vic_climatology_data(var.name=var.name,interval='1971-2000',
                                         clim=paste0(time,'_maximum'),
                                         read.dir=read.dir,
                                         gcm.list=gcm.list,ix=ix)

past.ens.avg <- calc(past.ens,mean)*1000

proj.ens <- extract_vic_climatology_data(var.name=var.name,interval='2041-2070',
                                       clim=paste0(time,'_maximum'),
                                       read.dir=read.dir,
                                       gcm.list=gcm.list,ix=ix)
proj.ens.avg <- calc(proj.ens,mean)*1000


if (1==0) {
###class.breaks <- c(0,25,50,100,150,200,250,300,400,500,1000,1500,2000)
class.breaks <- c(0,25,50,100,150,200,250,500,1000,1500,2000,2500,3000)

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/vic.ensemble.past.swe.peak.map.png')
plot.title <- 'VIC Peak SWE'

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')

make_van_whistler_plot(var.name='swe',plot.type='past',plot.title,past.ens.avg,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')
dev.off()

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/vic.ensemble.future.swe.peak.map.png')
plot.title <- 'VIC Peak SWE'

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')

make_van_whistler_plot(var.name='swe',plot.type='past',plot.title,proj.ens.avg,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')
dev.off()
}
##-----------------------------------------------------------------
##Anomalies

proj.ens.avg[proj.ens.avg > 3000] <- 3000
past.ens.avg[past.ens.avg > 3000] <- 3000
prct.ens.avg <- (proj.ens.avg - past.ens.avg)/past.ens.avg *100

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/vic.ensemble.percent.change.swe.peak.map.png')
plot.title <- 'VIC Peak SWE Anomalies (2041-2070 - 1971-2000)'
class.breaks <- seq(-80,0,10)

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')

make_van_whistler_plot(var.name='swe',plot.type='percent',plot.title,prct.ens.avg,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='%')
dev.off()
