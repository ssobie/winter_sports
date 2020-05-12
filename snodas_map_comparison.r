##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/ncc_map_support.r',chdir=T)
 

model <- 'ERA5'
reanalysis <- 'PNWNAmet'

if (1==1) {

  ##SNOW MODEL
  snow.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/calibrated_',model,'_',reanalysis,'_prism_tps/')
  ##snw.file <- paste0(snow.dir,'swe_BCCAQ2-PRISM_',model,'_',reanalysis,'_1945-2012.nc')
  snw.file <- paste0(snow.dir,'swe_BCCAQ2-PRISM_',model,'_',reanalysis,'_1980-2018.nc')
  snw.nc <- nc_open(snw.file)
  lon <- ncvar_get(snw.nc,'lon')
  lat <- ncvar_get(snw.nc,'lat')
  snow.time <- netcdf.calendar(snw.nc)

  ##SNODAS 
  snodas.dir <- paste0('/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/')
  snd.file <- paste0(snodas.dir,'swe_snodas_prism_grid_van_whistler_20100101-20181231.nc')
  snd.nc <- nc_open(snd.file)
  snodas.time <- netcdf.calendar(snd.nc)

  model.match <- format(snow.time,'%Y-%m-%d') %in% format(snodas.time,'%Y-%m-%d')
  snodas.match <- format(snodas.time,'%Y-%m-%d') %in% format(snow.time,'%Y-%m-%d')

  model.data <- array(NA,c(snw.nc$dim$lon$len,ncol=snw.nc$dim$lat$len,sum(model.match)))
 
  for (i in 1:snw.nc$dim$lon$len) {
     print(paste0(i,' of ',snw.nc$dim$lon$len))
     model.data[i,,] <- ncvar_get(snw.nc,'swe',start=c(i,1,1),count=c(1,-1,-1))[,model.match]*1000
  }
  
  snodas.raw <- ncvar_get(snd.nc,'swe')
  snodas.data <- snodas.raw[,,snodas.match]
  rm(snodas.raw)

  common.time <- snodas.time[snodas.match]

  data.diff <- apply(model.data - snodas.data,c(1,2),mean,na.rm=T)
  data.diff[data.diff > 1500] <- 1500
  data.diff[data.diff < -1500] <- -1500
  diff.raster <-  list(x=lon,y=lat,z=data.diff)
  save.dir <- '/storage/data/projects/rci/data/winter_sports/plots/data_files/'
  save(diff.raster,file=paste0(save.dir,model,'snodas.',model,'.',reanalysis,'.map.comparison.RData'))
}

##load(file=paste0(save.dir,model,'snodas.model.map.comparison.RData'))


class.breaks <- c(-1500,-1000,-750,-500,-250,-100,-50,0,50,100,250,500,750,1000,1500)
map.range <- range(data.diff,na.rm=T)
 ####get.class.breaks('swe',type='past',map.range,manual.breaks='')

plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/',model,'.prism.tps.swe.snodas.diff.2020.png') 
plot.title <- paste0(model,' Snow Model-VIC SWE Comparison')

png(file=plot.file,width=10,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='swe',plot.type='diff',plot.title,diff.raster,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')
dev.off()

browser()




##Need 
##ratio.raster for the success rate
##valid.raster for the observable days
##modis.snow.raster for the total snow days

save.dir <- '/storage/data/projects/rci/data/winter_sports/plots/data_files/'

##save(ratio.raster,file=paste0(save.dir,model,'.success.rate.2018.RData'))
##save(valid.raster,file=paste0(save.dir,'modis.observable.days.2018.RData'))
##save(modis.snow.raster,file=paste0(save.dir,'modis.snow.days.2018.RData'))
##save(model.snow.raster,file=paste0(save.dir,model,'.snow.days.2018.RData'))

nc_close(snd.nc)
nc_close(snw.nc)

##save.dir <- '/storage/data/projects/rci/data/winter_sports/plots/data_files/'
##load(paste0(save.dir,model,'.success.rate.2018.RData'))
##---------------------------------------------------------------------------------
##Comparison Plot
if (1==1)  {
plot.file <- paste('/storage/data/projects/rci/data/winter_sports/plots/',model,'.swe.snodas.diff.2020_manual2.png') 
plot.title <- paste0(model,' Snow Model-SNODAS SWE Comparison')

png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
map.range <- range(data.diff,na.rm=T)
leg.title <- 'mm'
class.breaks <- c(-10000,-1000,-500,-250,-100,-50,0,50,100,250,500,1000,100000)   ####get.class.breaks('swe',type='past',map.range,manual.breaks='')
map.class.breaks.labels <- get.class.break.labels(class.breaks,lesser.sign=TRUE,greater.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='swe',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)

vw.plot(diff.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE,bias=TRUE)
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
plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/modis.metro.van.comparison.data.2020.png'
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