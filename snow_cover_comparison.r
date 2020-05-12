##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/data/projects/rci/bcgov/moti/nrcan-precip_case_studies/code/moti.climdex.robjects.r',chdir=TRUE)
source('/storage/home/ssobie/code/hg/pievc/spatial.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/van_whistler_map.r')        

get.region.shape <- function(region,shape.dir) {
  region.shp <- readOGR(shape.dir, region, stringsAsFactors=F, verbose=F)
  return(region.shp)
}


##shape.dir <- paste0('/storage/data/projects/rci/data/assessments/metro_van/shapefiles/')
##region <- 'metro_van'
##region.shp <- spTransform(get.region.shape(region,shape.dir),CRS("+init=epsg:4326")) 

glacier.dir <- '/storage/data/gis/basedata/randolph_glacier_inventory/v32/02_rgi32_WesternCanadaUS'
glacier.name <- '02_rgi32_WesternCanadaUS'
glacier.shp <- spTransform(get.region.shape(glacier.name,glacier.dir),CRS("+init=epsg:4326"))

model <- 'ERA'
if (1==1) {
##---------------------------------------------------------------------------------
##MODIS 
modis.dir <- '/storage/data/projects/rci/data/winter_sports/MODIS_VAN_WHISTLER/'
snc.file <- paste0(modis.dir,'snc.modis.van_whistler.20010101-20181231.nc')
snc.nc <- nc_open(snc.file)
lon <- ncvar_get(snc.nc,'lon')
lat <- ncvar_get(snc.nc,'lat')
modis.time <- netcdf.calendar(snc.nc)
nc.grid <- get.netcdf.grid(snc.nc)
coordinates(nc.grid) <- c("lon", "lat")
modis.coords <- nc.grid@coords

##---------------------------------------------------------------------------------
##SNOW MODEL
snow.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/calibrated_manual/')
snw.file <- paste0(snow.dir,'snowdepth_BCCAQ2-PRISM_',model,'_19790101-20181031.nc')
snw.nc <- nc_open(snw.file)

lon <- ncvar_get(snw.nc,'lon')
lat <- ncvar_get(snw.nc,'lat')
snow.time <- netcdf.calendar(snw.nc)

nc.grid <- get.netcdf.grid(snw.nc)
coordinates(nc.grid) <- c("lon", "lat")
model.coords <- nc.grid@coords
date.match <- format(snow.time,'%Y-%m-%d') %in% format(modis.time,'%Y-%m-%d')
modis.match <- format(modis.time,'%Y-%m-%d') %in% format(snow.time,'%Y-%m-%d')
bnds <- range(which(date.match))
mbnds <- range(which(modis.match))


len <- snc.nc$dim$lon$len
cover.match <- matrix(NA,nrow=len,ncol=snc.nc$dim$lat$len)
cover.valid <- matrix(NA,nrow=len,ncol=snc.nc$dim$lat$len)

modis.snow <- matrix(NA,nrow=len,ncol=snc.nc$dim$lat$len)
model.snow <- matrix(NA,nrow=len,ncol=snc.nc$dim$lat$len)

for (i in 1:len) {
  print(i)
  snc.data <- ncvar_get(snc.nc,'snc',start=c(i,1,1),count=c(1,-1,-1))
  snc.modis <- snc.data[,mbnds[1]:mbnds[2]]
  snc.modis[snc.modis > 0 & snc.modis <=100] <- 1
  snc.modis[snc.modis > 100] <- NA
  snc.valid <- snc.modis
  ##snc.modis[snc.modis == 0]  <- NA

  snw.data <- ncvar_get(snw.nc,'snowdepth',start=c(i,1,bnds[1]),count=c(1,-1,diff(bnds)+1))
  snw.cover <- snw.data
  snw.cover[snw.cover>0.05] <- 1
  snw.cover[snw.cover<=0.05] <- 0
  
  snw.flagged <- snw.cover
  snw.flagged[is.na(snc.valid)] <- NA

  cover.diff <- snw.cover - snc.modis

  cover.match[i,] <- apply(cover.diff,1,function(x){sum(x==0,na.rm=T)})
  cover.valid[i,] <- apply(snc.valid,1,function(x){sum(!is.na(x))})

  modis.snow[i,] <- apply(snc.modis,1,sum,na.rm=T)
  model.snow[i,] <- apply(snw.flagged,1,sum,na.rm=T)

}

match.raster <- list(x=lon,y=lat,z=cover.match)
valid.raster <- list(x=lon,y=lat,z=cover.valid)
ratio.raster <- list(x=lon,y=lat,z=cover.match/cover.valid*100)
modis.snow.raster <- list(x=lon,y=lat,z=modis.snow)
model.snow.raster <- list(x=lon,y=lat,z=model.snow)

##Need 
##ratio.raster for the success rate
##valid.raster for the observable days
##modis.snow.raster for the total snow days

save.dir <- '/storage/data/projects/rci/data/winter_sports/plots/data_files/'

save(ratio.raster,file=paste0(save.dir,model,'.success.rate.2020.RData'))
save(valid.raster,file=paste0(save.dir,'modis.observable.days.2020.RData'))
save(modis.snow.raster,file=paste0(save.dir,'modis.snow.days.2020.RData'))
save(model.snow.raster,file=paste0(save.dir,model,'.snow.days.2020.RData'))

nc_close(snc.nc)
nc_close(snw.nc)

##browser()

}

save.dir <- '/storage/data/projects/rci/data/winter_sports/plots/data_files/'
load(paste0(save.dir,model,'.success.rate.2018.RData'))
##---------------------------------------------------------------------------------
##Comparison Plot
if (1==0)  {
plot.file <- paste('/storage/data/projects/rci/data/winter_sports/plots/',model,'.snc.success.rate.data.v2.png')
plot.title <- paste0(model,' Snow Model-MODIS Snow Cover Comparison')
map.range <- c(70,100)
leg.title <- 'Percent (%)'
class.breaks <- c(0,get.class.breaks('snc',type='past',map.range,manual.breaks=''))
map.class.breaks.labels <- get.class.break.labels(class.breaks,lesser.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='snd',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)

vw.plot(ratio.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE)
}

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
load(paste0(save.dir,'modis.snow.days.2018.RData'))
load(paste0(save.dir,'ERA.snow.days.2018.RData'))
model.snow <- model.snow.raster$z
modis.snow <- modis.snow.raster$z
lon <- model.snow.raster$x
lat <- model.snow.raster$y

plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/modis-era.percent.diff.2019.png'
plot.title <- '' ##MODIS Satellite Snow Cover Comparison (ERA - MODIS)'##paste0(model,'-MODIS Absolute Diff Sum')
leg.title <- 'Diff. (%)'
snow.prct <- (model.snow-modis.snow)/modis.snow*100
snow.prct[snow.prct > 100] <- 100
snow.prct[snow.prct < -100] <- -100
sp.raster <- list(x=lon,y=lat,z=snow.prct)
map.range <- range(snow.prct,na.rm=T)
##class.breaks <- seq(-500,500,100)
class.breaks <- c(-100,-75,-50,-20,-10,0,10,20,50,75,100) ##get.class.breaks('pr',type='anom',map.range,manual.breaks='') ##
map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=FALSE,lesser.sign=FALSE)
colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type='percent')
vw.plot(sp.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE)


##Success Rate for ERA
plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/modis-era.success.rate.2019.png'
plot.title <- '' ##paste0('ERA-Interim - MODIS All Days Success Rate')
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


}
##browser()
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
##Six panel figure
if (1==1) {
plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/modis.metro.van.comparison.data.2020_with_manual2.png'
png(plot.file,width=1800,height=1800)
par(mfrow=c(3,2))

##Observable days
plot.title <- 'MODIS Observable Days' 
load(paste0(save.dir,'modis.observable.days.2020.RData'))
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
load(paste0(save.dir,'modis.snow.days.2020.RData'))
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
plot.title <- paste0('ERA-Interim - MODIS All Days Success Rate')
load(paste0(save.dir,'ERA.success.rate.2020.RData'))
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
plot.title <- paste0('NCEP2 - MODIS All Days Success Rate')
load(paste0(save.dir,'NCEP2.success.rate.2020.RData'))
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


plot.title <- paste0('ERA-Interim - MODIS Snow Days Percent Difference')
load(paste0(save.dir,'ERA.snow.days.2020.RData'))
model.snow <- model.snow.raster$z
modis.snow <- modis.snow.raster$z
lon <- model.snow.raster$x
lat <- model.snow.raster$y
leg.title <- 'Diff. (days)'
snow.prct <- (model.snow-modis.snow)##/modis.snow*100
sp.raster <- list(x=lon,y=lat,z=snow.prct)
map.range <- range(snow.prct,na.rm=T)
class.breaks <- class.breaks <- c(-5000,seq(-500,500,50),5000) ##c(-5000,-100,-80,-60,-40,-20,-10,0,10,20,50,100,5000) ##c(-1000,-500,-250,-100,-50,-25,0,25,50,100,200,300,1000)  ##
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

plot.title <- paste0('NCEP2 - MODIS Snow Days Percent Difference')
load(paste0(save.dir,'NCEP2.snow.days.2020.RData'))
model.snow <- model.snow.raster$z
modis.snow <- modis.snow.raster$z
lon <- model.snow.raster$x
lat <- model.snow.raster$y

leg.title <- 'Diff. (days)'
snow.prct <- (model.snow-modis.snow) ##/modis.snow*100
sp.raster <- list(x=lon,y=lat,z=snow.prct)
map.range <- range(snow.prct,na.rm=T)
class.breaks <- c(-5000,seq(-500,500,50),5000) ##c(-5000,-100,-80,-60,-40,-20,-10,0,10,20,50,100,5000) ##c(-1000,-500,-250,-100,-50,-25,0,25,50,100,200,300,1000)  ##
map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=TRUE,lesser.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)
vw.plot(sp.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE)




dev.off()