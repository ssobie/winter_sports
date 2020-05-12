##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/ncc_map_support.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

library(raster)

model <- 'PNWNAmet'
reanalysis <- 'PNWNAmet'
type <- 'prism_tps'

##SNOW MODEL
snow.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/calibrated_',model,'_',reanalysis,'_',type,'/')
snw.file <- paste0(snow.dir,'swe_',tolower(model),'_con_swe.nc')
snw.nc <- nc_open(snw.file)
lon <- ncvar_get(snw.nc,'lon')
lat <- ncvar_get(snw.nc,'lat')
snow.time <- netcdf.calendar(snw.nc)

##VIC 
vic.dir <- paste0('/storage/data/projects/rci/data/winter_sports/')
vic.file <- paste0(vic.dir,'swe_day_VIC_BASE_historical_run1_19500101-20061231.nc')
vic.nc <- nc_open(vic.file)
vic.time <- netcdf.calendar(vic.nc)

model.match <- format(snow.time,'%Y-%m-%d') %in% format(vic.time,'%Y-%m-%d')
vic.match <- format(vic.time,'%Y-%m-%d') %in% format(snow.time,'%Y-%m-%d')

model.data <- ncvar_get(snw.nc,'swe')[,,model.match]*1000
vic.data <- ncvar_get(vic.nc,'swe')[,,vic.match]*1000

common.time <- vic.time[vic.match]
april1.ix <- grep('*-04-01',common.time)

data.diff <- apply(model.data - vic.data,c(1,2),mean,na.rm=T)
data.diff[data.diff > 500] <- 500
data.diff[data.diff < -500] <- -500

april.diff <- apply(model.data[,,april1.ix] - vic.data[,,april1.ix],c(1,2),mean,na.rm=T)

nc_close(vic.nc)
nc_close(snw.nc)

valid <- !is.na(data.diff)
vlen <- sum(valid)
vix <- which(valid,arr.ind=T)

cor.mat <- data.diff*0
sd.mat <- data.diff*0

for (i in 1:vlen) {
    print(i)
    ix <- as.numeric(vix[i,])
    cor.mat[ix[1],ix[2]] <- cor(model.data[ix[1],ix[2],],vic.data[ix[1],ix[2],])
    sd.model <- sd(model.data[ix[1],ix[2],])
    sd.vic <- sd(vic.data[ix[1],ix[2],])
    sd.mat[ix[1],ix[2]] <- sd.model/sd.vic
}


diff.raster <-  list(x=lon,y=lat,z=data.diff)
april.raster <-  list(x=lon,y=lat,z=april.diff)
cor.raster <-  list(x=lon,y=lat,z=cor.mat)
sd.raster <-  list(x=lon,y=lat,z=sd.mat)

##Need 
##ratio.raster for the success rate
##valid.raster for the observable days
##modis.snow.raster for the total snow days

save.dir <- '/storage/data/projects/rci/data/winter_sports/plots/data_files/'
##save(ratio.raster,file=paste0(save.dir,model,'.success.rate.2018.RData'))
##save(valid.raster,file=paste0(save.dir,'modis.observable.days.2018.RData'))
##save(modis.snow.raster,file=paste0(save.dir,'modis.snow.days.2018.RData'))
##save(model.snow.raster,file=paste0(save.dir,model,'.snow.days.2018.RData'))


##save.dir <- '/storage/data/projects/rci/data/winter_sports/plots/data_files/'
##load(paste0(save.dir,model,'.success.rate.2018.RData'))
##---------------------------------------------------------------------------------
##Average SWE Comparison Plot
if (1==1)  {

class.breaks <- c(-500,-300,-200,-100,-50,0,50,100,200,300,500)  
 ####get.class.breaks('swe',type='past',map.range,manual.breaks='')
plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/',model,'.',type,'.swe.vic.diff.2020.png')
plot.title <- paste0(model,' Snow Model-VIC SWE Comparison')

png(file=plot.file,width=10,height=6,units='in',res=600,pointsize=6,bg='white')
make_van_whistler_plot(var.name='swe',plot.type='diff',plot.title,diff.raster,plot.file,
                       class.breaks=class.breaks,mark=c(0,6,5,0),leg.title='mm')
dev.off()

browser()



##png(plot.file,width=1400,height=2700)
##par(mfrow=c(3,1))

map.range <- range(data.diff,na.rm=T)
leg.title <- 'mm'

map.class.breaks.labels <- get.class.break.labels(class.breaks,lesser.sign=TRUE,greater.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='swe',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)
plot.title <- paste0(model,' Snow Model-VIC Average SWE Comparison')
png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
vw.plot(diff.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE,bias=TRUE)
dev.off()



browser()




##SWE Corr
map.range <- range(cor.mat,na.rm=T)
leg.title <- 'Cor'
class.breaks <- c(0,0.25,0.5,0.6,0.7,0.8,0.9,1) ##get.class.breaks('swe',type='past',map.range,manual.breaks='')
map.class.breaks.labels <- get.class.break.labels(class.breaks,lesser.sign=FALSE,greater.sign=FALSE)
colour.ramp <- get.legend.colourbar(var.name='swe',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)
plot.title <- paste0(model,' Snow Model-VIC SWE Correlation')
vw.plot(cor.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE,bias=TRUE)


##SWE SD
map.range <- range(sd.mat,na.rm=T)
leg.title <- 'SD Ratio'
class.breaks <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,5,10000) ##get.class.breaks('swe',type='past',map.range,manual.breaks='')
map.class.breaks.labels <- get.class.break.labels(class.breaks,lesser.sign=FALSE,greater.sign=TRUE)
colour.ramp <- get.legend.colourbar(var.name='pr',map.range=map.range,
                                    my.bp=1,class.breaks=class.breaks,
                                    type)
plot.title <- paste0(model,' Snow Model-VIC Variability Ratio')
vw.plot(sd.raster,
        class.breaks,map.class.breaks.labels,colour.ramp,
        plot.file,plot.title,leg.title,
        glaciers=TRUE,bias=TRUE)

##dev.off()
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