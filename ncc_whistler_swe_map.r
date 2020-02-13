##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/ncc_whistler_mapping.r',chdir=T)

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
gcm.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/'
gcm.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')
time <- 'seasonal'
ix <- 1
var.name <- 'swe'
past.ens <- extract_climatology_data(var.name=var.name,interval='1971-2000',
                                       clim=paste0(time,'_mean'),
                                       read.dir=gcm.dir,
                                       gcm.list=gcm.list,ix=ix)
past.ens.avg <- calc(past.ens,mean)*1000

proj.ens <- extract_climatology_data(var.name=var.name,interval='2041-2070',
                                       clim=paste0(time,'_mean'),
                                       read.dir=gcm.dir,
                                       gcm.list=gcm.list,ix=ix)
proj.ens.avg <- calc(proj.ens,mean)*1000
proj.ens.avg[proj.ens.avg > 3000] <- 3000
past.ens.avg[past.ens.avg > 3000] <- 3000
prct.ens.avg <- (proj.ens.avg - past.ens.avg)/past.ens.avg *100

##-------------------------------------------------------------------------------------

snow.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/')


##---------------------------------------------------
##ERA Peak SWE

snow.file <- paste0(snow.dir,'swe_seasonal_mean_climatologies_BCCAQ2-PRISM_ERA_19790101-20181031.nc')
snow.seas <- brick(snow.file)*1000
snow.data <- subset(snow.seas,1)

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/whistler.swe.winter.mean.era.map.png')
plot.title <- 'Mean SWE'

##class.breaks <- c(0,25,50,100,150,200,250,500,1000,1500,2000,2500,3000)
class.breaks <- c(0,25,50,100,150,200,250,300,400,500,1000,1500,2000)

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_whistler_plot(var.name='swe',plot.type='past',plot.title,snow.data,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')


dev.off()


##---------------------------------------------------
##Model Past Peak SWE

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/whistler.swe.winter.mean.ENS.1971-2000.png')
plot.title <- 'Peak SWE (1971-2000)'

##class.breaks <- c(0,25,50,100,150,200,250,500,1000,1500,2000,2500,3000)
class.breaks <- c(0,25,50,100,150,200,250,300,400,500,1000,1500,2000)

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_whistler_plot(var.name='swe',plot.type='past',plot.title,past.ens.avg,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')

dev.off()


##---------------------------------------------------
##Model Future Peak SWE

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/whistler.swe.winter.mean.ENS.2041-2070.png')
plot.title <- 'Peak SWE (2041-2070)'

##class.breaks <- c(0,25,50,100,150,200,250,500,1000,1500,2000,2500,3000)
class.breaks <- c(0,25,50,100,150,200,250,300,400,500,1000,1500,2000)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_whistler_plot(var.name='swe',plot.type='past',plot.title,proj.ens.avg,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')


dev.off()


##---------------------------------------------------
##Model Percent Peak SWE

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/whistler.percent.swe.winter.mean.ENS.2041-2070.1971-2000.png')
plot.title <- 'Percent Change in Peak SWE (2041-2070 - 1971-2000)'
##class.breaks <- c(0,25,50,75,100,150,200,250,300,350,400,500)
class.breaks <- NULL ##c(0,25,50,100,150,200,250,500,1000,1500,2000,2500,3000)

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_whistler_plot(var.name='swe',plot.type='percent',plot.title,prct.ens.avg,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='%')

dev.off()
