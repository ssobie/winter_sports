##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)

source('/storage/home/ssobie/code/repos/winter_sports/test.snow.model.r',chdir=T)

get.coordinates <- function(site) {

  coordinates <- list(callaghan=c(-123.1036,50.1383278,1009),
                      orchid_lake=c(-123.0519638,49.53678,1178),
                      palisade_lake=c(-123.0321944,49.454433,898),
                      grouse_mountain=c(-123.0774472,49.383655,1126),
                      dog_mountain=c(-122.96255,49.37251944,1007),
                      dickson_lake=c(-122.06984166,49.3168194,1147),
                      stave_lake=c(-122.315805,49.58030277,1211),
                      nahatlatch=c(-122.059261,49.825866,1530),
                      wahleach=c(-121.57945,49.2298694,1395),
                      klesilkwa=c(-121.3086527,49.129438,610),
                      lightning_lake=c(-120.850205,49.044788,1254),
                      brookmere=c(-120.87397,49.815027,994),
                      shovelnose_mountain=c(-120.864175,49.8546305,1456),
                      hamilton_hill=c(-120.7955805,49.4988027,1477),
                      spuzzum_creek=c(-121.686,49.674,1197),
                      chilliwack_river=c(-121.71667,49.0333,1600),
                      upper_squamish=c(-123.4333,50.1500,1340),
                      wahleach_lake=c(-121.5833,49.2333,1400),
                      tenquille_lake=c(-122.9333,50.5333,1680))

  rv <- coordinates[[site]]
  return(rv)
}

hyper.snow <- function(pr,tasmax,tasmin,coeffs) {
        
        tas <- (tasmax+tasmin)/2                # degrees C
        frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)
        sample <- runif(length(tas),min=0,max=100)
        test <- sample > frac        
        snow.type <- rep(TRUE,length(tas))
        snow.type[test] <- FALSE

        NewSnowWatEq <- pr*0.001
        NewSnowWatEq[!snow.type] <- 0
        R_m <- pr*0.001
        R_m[snow.type] <- 0
        rv <- list(snow=NewSnowWatEq,
                   rain=R_m)
        return(rv)
}


##Slope and Aspect Values
as.dir <- '/storage/data/projects/rci/data/prism/'
slopes.nc <- nc_open(paste0(as.dir,'prism_slopes.nc'))
bc.slopes <- ncvar_get(slopes.nc,'Band1')/90*pi/2
bc.lon <- ncvar_get(slopes.nc,'lon')
bc.lat <- ncvar_get(slopes.nc,'lat')
nc_close(slopes.nc)

aspects.nc <- nc_open(paste0(as.dir,'prism_aspects.nc')) 
bc.aspects <- ncvar_get(aspects.nc,'Band1')/360*2*pi
nc_close(aspects.nc)



sites <- c('callaghan',
           'orchid_lake',
           'palisade_lake',
           'grouse_mountain',
           'dog_mountain',
           'stave_lake',
           'nahatlatch',
           'wahleach',
           'klesilkwa',
           'lightning_lake',
           'brookmere',
           'shovelnose_mountain',
           'hamilton_hill')
sites <- 'shovelnose_mountain'
model <- 'NCEP2'

course.site.swe <- vector(mode='list',length=length(sites))
model.site.swe <- vector(mode='list',length=length(sites))

##sites <- c('callaghan','orchid_lake')
slen <- 10
swe.sims <- matrix(0,nrow=slen,ncol=13696)
snow.sims <- matrix(0,nrow=slen,ncol=13696)

##Loop over sites

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)
    ##Reanalysis 800m data
    clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sites/',site,'_',model,'_800m_data.csv')
    clim.data <- read.csv(clim.file,header=T,as.is=T)
    pr.data <- clim.data$Pr
    tasmax.data <- clim.data$Tasmax
    tasmin.data <- clim.data$Tasmin
    tas.data <- clim.data$Tas
    dates <- clim.data$Dates

    coords <- get.coordinates(site)
    lat.bnds <- coords[2]
    elev <- coords[3]
    site.slope <- bc.slopes[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]    
    site.aspect <- bc.aspects[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]    

    ##Snow Course Data
    course.file <- paste('/storage/data/projects/rci/data/assessments/snow_model/',site,'_snow_course.csv',sep='')
    course.data <- read.csv(course.file,header=T,as.is=T)
    course.dates <- format(as.Date(course.data[,1]),'%Y-%m-%d')
    course.swe <- course.data[,3] ##mm
    course.pack <- course.data[,2] ##cm
    course.dense <-  course.data[,4]
    sb <- 1:length(course.dates) ##190:210
    date.subset <- format(as.Date(dates),'%Y-%m-%d') %in% course.dates
    course.subset <- course.dates %in% format(as.Date(dates),'%Y-%m-%d')

##    course.site.swe[[i]] <- course.swe[course.subset]
##    model.site.swe[[i]]  <- results$swe[date.subset]*1000

    ##coeffs <- list(a=-48.2372,b=0.7449,c=1.0919,d=1.0209)
    ##coeffs <- list(a=-49.49,b=0.4,c=4.5,d=1.0209)
##    coeffs <- list(a=-49.49,b=0.4128,c=2.6545,d=1.0209)
    coeffs <- list(a=-49.49,b=0.5628,c=1.5,d=1.0209)

    for (k in 1:slen) {
    print(k)
    test.snow <- hyper.snow(pr.data,tasmax.data,tasmin.data,coeffs)
    results <- snow.melt(Date=dates, precip_mm=pr.data, Tmax_C=tasmax.data, Tmin_C=tasmin.data,Snow=test.snow,
                         lat_deg=lat.bnds, slope=site.slope, aspect=site.aspect, tempHt=1, windHt=2, groundAlbedo=0.25,
                         SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)
    swe.sims[k,] <- results$swe
    snow.sims[k,] <- results$snowdepth
    }                                      

if (1==0) {
    ##Snow Pillow Data
    ##pillow.file <- paste('/storage/data/projects/rci/data/assessments/snow_model/snow_pillow/',site,'_asp.csv',sep='')
    ##pillow.data <- read.csv(pillow.file,header=T,as.is=T)
    ##pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
    ##pillow.tasmax <- pillow.data[,3]
    ##pillow.tasmin <- pillow.data[,5]
    ##pillow.tas <- (pillow.tasmax + pillow.tasmin)/2
    ##pillow.precip <- pillow.data[,7]##mm
    ##pillow.swe <- pillow.data[,11] ##mm
    ##pillow.pack <- pillow.data[,13] ##cm
    ##sb <- 2800:3200
    ##sb <- 1:length(pillow.dates)

##    asp.model <- snow.melt(Date=pillow.dates[sb], precip_mm=pillow.precip[sb], Tmax_C=pillow.tasmax[sb], Tmin_C=pillow.tasmin[sb], 
##                         lat_deg=lat.bnds, slope=site.slope, aspect=site.aspect, tempHt=1, windHt=2, groundAlbedo=0.25,
##                         SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)


    date.subset <- format(as.Date(dates),'%Y-%m-%d') %in% pillow.dates[sb]

    plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
    png(file=paste0(plot.dir,model,'.',site,'.snow.pillow.comparison.png'),width=800,height=900)
    par(mfrow=c(3,1))    

    ##plot(as.Date(dates[date.subset]),round(tasmin.data[date.subset],1),type='l',lwd=3,col='red',main='TASMIN (C)',cex.axis=1.5)
    ##lines(as.Date(pillow.dates),pillow.tasmin,lwd=3,col='orange')

    plot(as.Date(dates[date.subset]),round(tas.data[date.subset],1),type='l',lwd=3,col='red',main=paste0(site,' ',model,'\nTAS (C)'),cex.axis=1.5,
         xlab='Date',ylab='Average Temperature (C)', cex.lab=1.5)
    lines(as.Date(pillow.dates),pillow.tas,lwd=3,col='orange')
    abline(h=0)
##    plot(as.Date(dates[date.subset]),round(tas.diff,1),type='l',lwd=3,col='red',main='TAS (C)',cex.axis=1.5)

##    plot(as.Date(dates[date.subset]),round(tas.data[date.subset],1),type='l',lwd=3,col='red',main='DS TAS - Pillow TAS (C)',cex.axis=1.5)
##    lines(as.Date(pillow.dates),pillow.tas,lwd=3,col='orange')
##    abline(h=0)
##    plot(as.Date(dates[date.subset]),round(tas.diff,1),type='l',lwd=3,col='red',main='TAS (C)',cex.axis=1.5)

    plot(as.Date(dates[date.subset]),(pr.data[date.subset]),type='l',lwd=3,col='blue',main='Precip (mm)',cex.axis=1.5,
         xlab='Date',ylab='Precipitation (mm)', cex.lab=1.5)
    lines(as.Date(pillow.dates)[sb],(pillow.precip[sb]),lwd=3,col='green')
    abline(h=0)


##    plot(as.Date(dates[date.subset]),round(results$snowfall[date.subset]*100,2),type='l',lwd=3,col='blue',main='Snowfall (cm)',cex.axis=1.5)
##    lines(as.Date(pillow.dates),pillow.precip,lwd=3,col='green')
                   

    plot(as.Date(dates)[date.subset],results$snowdepth[date.subset]*100,type='l',lwd=3,col='blue',ylim=c(0,900),cex.axis=1.5,
             xlab='Date',ylab='Snowpack (cm)', cex.lab=1.5)
    points(as.Date(pillow.dates),pillow.pack,cex=1,col='green',pch=16)
    abline(h=0)
##    lines(as.Date(pillow.dates[sb]),asp.model$snowdepth*100,lwd=3,col='yellow')


##    pillow.mon.fac <- as.factor(format(as.Date(pillow.dates),'%Y-%m'))
##    pillow.pr.mon <- tapply(pillow.precip,pillow.mon.fac,sum,na.rm=T)
##    pillow.sd.mon <- tapply(pillow.pack,pillow.mon.fac,mean,na.rm=T)

##    model.mon.fac <- as.factor(format(as.Date(dates[date.subset]),'%Y-%m'))
##    model.pr.mon <- tapply(pr.data[date.subset],model.mon.fac,sum,na.rm=T)
##    model.sd.mon <- tapply(results$snowdepth[date.subset],model.mon.fac,sum,na.rm=T)
    dev.off()

}

    ds.mon.fac <- as.factor(format(as.Date(dates),'%Y-%m'))
    ds.mon.dates <- paste0(levels(ds.mon.fac),'-01')
    pr.mon.tot <- tapply(pr.data,ds.mon.fac,sum,na.rm=T)

##Snow Course Comparison

    plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/course_comparison/'
    png(file=paste0(plot.dir,model,'.',site,'.snow.course.comparison.png'),width=1000,height=600)
    par(mfrow=c(2,2))    

    plot(as.Date(ds.mon.dates)[sb],pr.mon.tot[sb],cex=1,col='blue',lwd=3,type='l',
           cex.axis=1.5,main='Monthly Precipitation (mm)',
           xlab='Date',ylab='Precip (mm)', cex.lab=1.5)
    abline(h=0)

    plot(as.Date(course.dates)[sb],course.pack[sb],cex=1,col='green',pch=16,ylim=c(0,max(results$snowdepth*100)),cex.axis=1.5,main='Snow Pack',
             xlab='Date',ylab='Snowpack (cm)', cex.lab=1.5)
    apply(snow.sims,1,function(x,y){lines(y,x*100,col='lightblue',lwd=2)},as.Date(dates))
    lines(as.Date(dates),results$snowdepth*100,lwd=3,col='blue')
    points(as.Date(course.dates),course.pack,cex=1,col='green',pch=16)
    abline(h=0)

    plot(as.Date(course.dates)[sb],course.swe[sb],cex=1,col='green',pch=16,ylim=c(0,max(results$swe*1000)),cex.axis=1.5,main='SWE',
             xlab='Date',ylab='SWE (mm)', cex.lab=1.5)
    apply(swe.sims,1,function(x,y){lines(y,x*1000,col='lightblue',lwd=2)},as.Date(dates))
    lines(as.Date(dates),results$swe*1000,lwd=3,col='blue')
    points(as.Date(course.dates),course.swe,cex=1,col='green',pch=16)
    abline(h=0)


    plot(as.Date(course.dates)[sb],course.dense[sb],cex=1,col='green',pch=16,ylim=c(0,70),cex.axis=1.5,main='Snow Density',
             xlab='Date',ylab='Snow Density (%)', cex.lab=1.5)
    lines(as.Date(dates),results$snowdense,lwd=3,col='blue')
    points(as.Date(course.dates),course.dense,cex=1,col='green',pch=16)
    abline(h=0)
    dev.off()
    ##taylor.diagram(course.pack[course.subset],results$snowdepth[date.subset]*1000,add=TRUE)
        
}        





##Compare ERA, NCEP2 snow depth against course depth

##Compare ERA, NCEP2 SWE against course SWE

##Compare ERA, NCEP2 snow density against course density

##Compare ERA, NCEP2 snow cover against MODIS snow cover

##Show group plots of all sites
