##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)

source('/storage/home/ssobie/code/repos/winter_sports/test.snow.model.r',chdir=T)

get.coordinates <- function(site) {

  coordinates <- list(spuzzum_creek=c(-121.686,49.674,1197),
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



model <- 'ERA'


sites <- c('upper_squamish') ##,'spuzzum_creek','tenquille_lake')
site.name <- 'Upper Squamish'
slen <- 10

pillow.site.swe <- vector(mode='list',length=length(sites))
model.site.swe <- vector(mode='list',length=length(sites))

swe.sims <- matrix(0,nrow=slen,ncol=13696)
snow.sims <- matrix(0,nrow=slen,ncol=13696)

##Loop over sites

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)
    ##Reanalysis 800m data

    ##clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sites/sp_testing/',site,'_',model,'_800m_data_196_55.csv')
    clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sites/',site,'_',model,'_800m_data.csv')
    clim.data <- read.csv(clim.file,header=T,as.is=T)
    pr.data <- clim.data$Pr
    tasmax.data <- clim.data$Tasmax
    tasmin.data <- clim.data$Tasmin
    tas.data <- clim.data$Tas
    dates <- clim.data$Dates
if (1==0) {
    coords <- get.coordinates(site)
    lat.bnds <- coords[2]
    elev <- coords[3]
    site.slope <- bc.slopes[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]    
    site.aspect <- bc.aspects[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]    

    ##coeffs <- list(a=-49.49,b=0.4128,c=2.6545,d=1.0209)
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
    }

    swe.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sims/',site,'.',tolower(model),'.swe.1001.csv')
    swe.sims <- read.csv(swe.file,header=T,as.is=T)
    snow.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sims/',site,'.',tolower(model),'.snow.1001.csv')
    snow.sims <- read.csv(snow.file,header=T,as.is=T)


    winter.flag <- grep('(*-12-*|*-01-*|*-02-*)',dates)
    swe.winter <- swe.sims[winter.flag,]
    swe.winter.pos <- swe.winter[swe.winter>0]

    pack.winter <- snow.sims[winter.flag,]
    pack.winter.pos <- pack.winter[pack.winter>0]

    ##Snow Pillow Data
    pillow.file <- paste('/storage/data/projects/rci/data/assessments/snow_model/snow_pillow/',site,'_asp.csv',sep='')
    pillow.data <- read.csv(pillow.file,header=T,as.is=T)
    pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
    pillow.tasmax <- pillow.data[,3]
    pillow.tasmin <- pillow.data[,5]
    pillow.tas <- (pillow.tasmax + pillow.tasmin)/2
    pillow.precip <- pillow.data[,7]##mm
    pillow.swe <- pillow.data[,11] ##mm
    pillow.pack <- pillow.data[,13] ##cm
    
    sb <- 1850:7150
    
    date.subset <- format(as.Date(dates),'%Y-%m-%d') %in% pillow.dates[sb]
    pillow.subset <- pillow.dates %in% format(as.Date(dates[date.subset]),'%Y-%m-%d')

    ##Snow season start and end dates
    dates.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/snow_seasons/',site,'_season_dates.csv')
    season.dates <- read.csv(dates.file,header=T,as.is=T)
    pillow.starts <- as.Date(season.dates$Start)
    pillow.ends <- as.Date(season.dates$End)
    pillow.lengths <- as.numeric(season.dates$Length)


    pillow.winter.flag <- grep('(*-12-*|*-01-*|*-02-*)',pillow.dates)
    pillow.swe.winter <- pillow.swe[pillow.winter.flag]    
    p.swe.winter.pos <- pillow.swe.winter[pillow.swe.winter>0]

    pillow.pack.winter <- pillow.pack[pillow.winter.flag]    
    p.pack.winter.pos <- pillow.pack.winter[pillow.pack.winter>0]


    breaks <- seq(0,3000,100)
    m.breaks <- seq(0,2900,100)    
    p.breaks <- seq(50,2950,100)
    m.hist <- hist(swe.winter.pos,breaks=breaks,plot=F)
    p.hist <- hist(p.swe.winter.pos,breaks=breaks,plot=F)

    mp.breaks <- seq(0,750,50)
    m.p.breaks <- seq(0,700,50)
    p.p.breaks <- seq(25,725,50)
    m.p.hist <- hist(pack.winter.pos,breaks=mp.breaks,plot=F)
    p.p.hist <- hist(p.pack.winter.pos,breaks=mp.breaks,plot=F)

    ##plot(as.Date(dates[date.subset]),round(tasmin.data[date.subset],1),type='l',lwd=3,col='red',main='TASMIN (C)',cex.axis=1.5)
    ##lines(as.Date(pillow.dates),pillow.tasmin,lwd=3,col='orange')

        
    plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
    png(file=paste0(plot.dir,model,'.',site,'.snow.pillow.histogram.png'),width=1200,height=900)
    par(mar=c(6,5,5,3))      
    par(mfrow=c(2,2))    

    ##Winter SWE Histogram
    plot(x=m.breaks+20,y=m.hist$density*100,type='h',lwd=10,xlab='SWE (mm)',ylab='Fraction (%)',main='Snow Water Equivalent',lend=2,
         cex.lab=1.75,cex.main=2,cex.axis=1.75,yaxs='i',col='blue',ylim=c(0,0.18))
    lines(x=p.breaks+15,y=p.hist$density*100,type='h',lwd=10,col='black',lend=2)    

    ##Winter Snowpack Histogram
    plot(x=m.p.breaks+12.5,y=m.p.hist$density*50,type='h',lwd=20,ylim=c(0,0.30),yaxs='i',col='blue',
         xlab='Snowdepth (cm)',ylab='Fraction (%)',main='Snowpack',lend=2,cex.lab=1.75,cex.main=2,cex.axis=1.75)
    lines(x=p.p.breaks+10,y=p.p.hist$density*50,type='h',lwd=20,col='black',lend=2)    
    legend('topright',leg=c('Obs.','Model'),col=c('black','blue'),pch=15,cex=2)
    ##SWE
    plot(as.Date(dates[date.subset]),round(swe.sims[date.subset,1],2),xlab='Date',ylab='SWE (mm)',yaxs='i',
         type='l',lwd=3,col='white',main='Snow Water Equivalent',cex.axis=1.75,cex.lab=1.75,cex.main=2,ylim=c(0,3000))
    apply(swe.sims[date.subset,],2,function(x,y){lines(y,x,col='lightblue',lwd=2)},as.Date(dates[date.subset]))
    points(as.Date(pillow.dates),pillow.swe,pch=16,col='black')
    lines(as.Date(dates[date.subset]),apply(swe.sims,1,mean)[date.subset],col='blue',lwd=3)
    ##abline(v=as.Date(dates[date.subset][227]))
                  
    ##Snowpack                   
    plot(as.Date(dates)[date.subset],snow.sims[date.subset,1],type='l',lwd=3,col='white',ylim=c(0,800),cex.axis=1.5,yaxs='i',
             xlab='Date',ylab='Snowpack (cm)',main='Snowpack', cex.lab=1.75,cex.axis=1.75,cex.main=2)
    apply(snow.sims[date.subset,],2,function(x,y){lines(y,x,col='lightblue',lwd=2)},as.Date(dates[date.subset]))
    points(as.Date(pillow.dates),pillow.pack,cex=1,col='black',pch=16)
    lines(as.Date(dates[date.subset]),apply(snow.sims,1,mean)[date.subset],col='blue',lwd=3)
    abline(h=0)
    legend('topright',leg=c('Obs.','Model'),col=c('black','blue'),pch=15,cex=2)
    ##abline(v=as.Date(dates[date.subset][227]))
    dev.off()
browser()
    ##Temperature
    
##    plot(as.Date(dates[date.subset]),round(tas.data[date.subset],1),type='l',lwd=3,col='red',
##         main=paste0(site,' ',model,'\nTAS (C)'),cex.axis=1.5,
##         xlab='Date',ylab='Average Temperature (C)', cex.lab=1.5)
    plot(as.Date(dates[date.subset]),round(tas.data[date.subset] - pillow.tas[pillow.subset],1),type='l',lwd=3,col='red',
         main=paste0(site,' ',model,'\nTAS (C)'),cex.axis=1.5,
         xlab='Date',ylab='Average Temperature (C)', cex.lab=1.5)
##    lines(as.Date(pillow.dates),pillow.tas,lwd=3,col='orange')
##    lines(as.Date(dates[date.subset]),tas.data[date.subset]-pillow.tas[pillow.subset],lwd=3,col='green')
    abline(h=0)
    abline(v=as.Date(dates[date.subset][227]))

    ##Precipitation
    
##    plot(as.Date(dates[date.subset]),(pr.data[date.subset]),type='l',lwd=3,col='blue',main='Precip (mm)',cex.axis=1.5,
##         xlab='Date',ylab='Precipitation (mm)', cex.lab=1.5) 
##    lines(as.Date(pillow.dates)[sb],(pillow.precip[sb]),lwd=3,col='green')
##    points(as.Date(dates[date.subset]),pillow.precip[pillow.subset]-pr.data[date.subset],pch='-',cex=2,col='red')    
##    abline(h=0)
##    abline(v=as.Date(dates[date.subset][227]))
##    pr.diff <- pillow.precip[pillow.subset]-pr.data[date.subset]        
##    pr.diff[is.na(pr.diff)] <- 0

    
    cumulative.pr <- pillow.precip*0

    for (y in 1:nrow(season.dates)) {
  
      csp.subset <- c(grep(pillow.starts[y],pillow.dates):grep(pillow.ends[y],pillow.dates))
      csm.subset <- c(grep(pillow.starts[y],dates):grep(pillow.ends[y],dates))
    
      csp.precip <- pillow.precip[csp.subset]
      csm.precip <- pr.data[csm.subset]
    
      csp.flag <- is.na(csp.precip)
      csp.precip[csp.flag] <- 0
      csm.precip[csp.flag] <- 0
      cp.precip.diff <- csm.precip - csp.precip
      tp.precip <- cumsum(cp.precip.diff)
    
      cumulative.pr[csp.subset] <- tp.precip 
    
    }

    plot(as.Date(pillow.dates[sb]),(cumulative.pr[sb]),cex=2,col='blue',main='Precip Diff (mm)',cex.axis=1.5,
         xlab='Date',ylab='Precipitation (mm)', cex.lab=1.5,pch='-') 
    ##lines(as.Date(pillow.dates)[sb],(pillow.precip[sb]),lwd=3,col='green')
    ##points(as.Date(dates[date.subset]),pillow.precip[pillow.subset]-pr.data[date.subset],pch='-',cex=2,col='red')    
    abline(h=0)
    abline(v=as.Date(dates[date.subset][227]))
    pr.diff <- pillow.precip[pillow.subset]-pr.data[date.subset]        
    pr.diff[is.na(pr.diff)] <- 0
    

}