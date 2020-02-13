##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/ncc_map_support.r',chdir=T)

library(raster)


extract_ens_climatology_data <- function(var.name,type,interval,clim,
                                         read.dir) {

   ens.dir <- paste0(read.dir,type,'/ENSEMBLE/')
   var.files <- list.files(path=ens.dir,pattern=paste0(var.name,'_',clim,'_climatology'))
   ens.file <- var.files[grep(interval,var.files)]
   ens.raster <- subset(brick(paste0(ens.dir,ens.file)),1)
   return(ens.raster)
}


if (1==0) {

prism.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/ERA/')

##---------------------------------------------------
#Winter Minimum Temperature

prism.file <- paste0(prism.dir,'tasmin_seasonal_climatologies_ERA_1980-2018.nc')
prism.seas <- brick(prism.file)
prism.winter <- subset(prism.seas,1)

prism.data <- prism.winter
plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/prism.winter.tasmin.map.png')
plot.title <- 'Winter Minimum Temperature'
class.breaks <- c(-15,-10,-8,-6,-4,-2,-1,0,1,2,4)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='tasmin',plot.type='past',plot.title,prism.data,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='\u00B0C')
dev.off()

##---------------------------------------------------
#Winter Total Precipitation

prism.file <- paste0(prism.dir,'pr_seasonal_total_climatologies_ERA_1980-2018.nc')
prism.seas <- brick(prism.file)
prism.winter <- subset(prism.seas,1)

prism.data <- prism.winter
plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/prism.winter.pr.map.png')
plot.title <- 'Winter Average Precipitation'
class.breaks <- c(0,250,500,600,700,800,900,1000,1250,1500,1750,2250)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='pr',plot.type='past',plot.title,prism.data,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')
dev.off()

##---------------------------------------------------
}



##-----------------------------------------------------
##Future projections
past.int <- '1971-2000'
proj.int <- '2041-2070'
##-------------------------
read.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/bc/rcp85/'

prw.past.ens <- extract_ens_climatology_data(var.name='pr',type='seasonal',interval=past.int,
                                     clim='winter_average',read.dir=read.dir)
prw.proj.ens <- extract_ens_climatology_data(var.name='pr',type='seasonal',interval=proj.int,
                                     clim='winter_average',read.dir=read.dir)
prw.prct.ens <- (prw.proj.ens - prw.past.ens)/prw.past.ens * 100 
prw.anom.ens <- prw.proj.ens - prw.past.ens



plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/anomaly.change.winter.pr.prism.map.png')
plot.title <- 'Winter Average Precipitation (2041-2070 - 1971-2000)'
class.breaks <- c(0,10,20,30,40,50,60,70,80)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='pr',plot.type='past',plot.title,prw.anom.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')
dev.off()


plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/percent.change.winter.pr.prism.map.png')
plot.title <- 'Winter Average Precipitation (2041-2070 - 1971-2000)'
class.breaks <- seq(0,14,2)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='pr',plot.type='past',plot.title,prw.prct.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='%')
dev.off()

browser()


prs.past.ens <- extract_ens_climatology_data(var.name='pr',type='seasonal',interval=past.int,
                                     clim='summer_average',read.dir=read.dir)
prs.proj.ens <- extract_ens_climatology_data(var.name='pr',type='seasonal',interval=proj.int,
                                     clim='summer_average',read.dir=read.dir)
prs.prct.ens <- (prs.proj.ens - prs.past.ens)/prs.past.ens * 100 
prs.anom.ens <- prs.proj.ens - prs.past.ens


plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/anomaly.change.summer.pr.prism.map.png')
plot.title <- 'Summer Average Precipitation (2041-2070 - 1971-2000)'
class.breaks <- seq(-120,0,20)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='pr',plot.type='past',plot.title,prs.anom.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')
dev.off()


plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/percent.change.summer.pr.prism.map.png')
plot.title <- 'Summer Average Precipitation (2041-2070 - 1971-2000)'
class.breaks <- seq(-25,-5,by=5)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='pr',plot.type='past',plot.title,prs.prct.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='%')
dev.off()

##-------------------------
##Spring
read.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/bc/rcp85/'
past.int <- '1971-2000'
proj.int <- '2041-2070'

prp.past.ens <- extract_ens_climatology_data(var.name='pr',type='seasonal',interval=past.int,
                                     clim='spring_average',read.dir=read.dir)
prp.proj.ens <- extract_ens_climatology_data(var.name='pr',type='seasonal',interval=proj.int,
                                     clim='spring_average',read.dir=read.dir)
prp.prct.ens <- (prp.proj.ens - prp.past.ens)/prp.past.ens * 100 
prp.anom.ens <- prp.proj.ens - prp.past.ens

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/past.spring.pr.prism.map.png')
plot.title <- 'Spring Average Precipitation 1971-2000'
##class.breaks <- NULL #seq(-120,0,20)
class.breaks <- c(0,250,500,600,700,800,900,1000,1250,1500,1750,2250)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='pr',plot.type='past',plot.title,prp.past.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')
dev.off()


plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/anomaly.change.spring.pr.prism.map.png')
plot.title <- 'Spring Average Precipitation (2041-2070 - 1971-2000)'
class.breaks <- NULL #seq(-120,0,20)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='pr',plot.type='past',plot.title,prp.anom.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')
dev.off()


plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/percent.change.spring.pr.prism.map.png')
plot.title <- 'Spring Average Precipitation (2041-2070 - 1971-2000)'
class.breaks <- seq(0,14,by=2)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='pr',plot.type='past',plot.title,prp.prct.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='%')
dev.off()

browser()


##-------------------------
if (1==0) {
tnw.past.ens <- extract_ens_climatology_data(var.name='tasmin',type='seasonal',interval=past.int,
                                     clim='winter_average',read.dir=read.dir)
tnw.proj.ens <- extract_ens_climatology_data(var.name='tasmin',type='seasonal',interval=proj.int,
                                     clim='winter_average',read.dir=read.dir)
tnw.anom.ens <- tnw.proj.ens - tnw.past.ens
 

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/percent.change.winter.tasmin.prism.map.png')
plot.title <- 'Winter Average Temperature (2041-2070 - 1971-2000)'
class.breaks <- seq(1.8,3.6,0.2)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='tasmin',plot.type='past',plot.title,tnw.anom.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='\u00B0C')
dev.off()


plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/past.winter.tasmin.prism.map.png')
plot.title <- 'Winter Average Temperature (1971-2000)'
class.breaks <- c(-15,-12.5,-10,-7.5,-5,-3,-2,-1,0,1,2,3,5)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='tasmin',plot.type='past',plot.title,tnw.past.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='\u00B0C')
dev.off()


plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/future.winter.tasmin.prism.map.png')
plot.title <- 'Winter Average Temperature (2041-2070)'
class.breaks <- c(-15,-12.5,-10,-7.5,-5,-3,-2,-1,0,1,2,3,5)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='tasmin',plot.type='past',plot.title,tnw.proj.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='\u00B0C')
dev.off()

}
##-------------------------------
##Spring Minimum Temperature

tnp.past.ens <- extract_ens_climatology_data(var.name='tasmin',type='seasonal',interval=past.int,
                                     clim='spring_average',read.dir=read.dir)
tnp.proj.ens <- extract_ens_climatology_data(var.name='tasmin',type='seasonal',interval=proj.int,
                                     clim='spring_average',read.dir=read.dir)
tnp.anom.ens <- tnp.proj.ens - tnp.past.ens
 

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/anomaly.change.spring.tasmin.prism.map.png')
plot.title <- 'Spring Average Temperature (2041-2070 - 1971-2000)'
class.breaks <- NULL
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='tasmin',plot.type='past',plot.title,tnp.anom.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='\u00B0C')
dev.off()


plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/past.spring.tasmin.prism.map.png')
plot.title <- 'Spring Average Temperature (1971-2000)'
class.breaks <- seq(-20,10,2.5)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='tasmin',plot.type='past',plot.title,tnp.past.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='\u00B0C')
dev.off()


plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/future.spring.tasmin.prism.map.png')
plot.title <- 'Spring Average Temperature (2041-2070)'
class.breaks <- seq(-20,10,2.5)
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='tasmin',plot.type='past',plot.title,tnp.proj.ens,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='\u00B0C')
dev.off()
