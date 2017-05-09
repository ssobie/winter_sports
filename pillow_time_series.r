##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)

source('/storage/home/ssobie/code/repos/winter_sports/test.snow.model.r',chdir=T)

get.coordinates <- function(site) {

  coordinates <- list(spuzzum_creek=c(-121.686,49.674,1197),
                      chilliwack_river=c(-121.71667,49.0333,1600),
                      upper_squamish=c(-123.4333,50.1500,1340),
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



model <- 'NCEP2'
type <- 'SWE'


sites <- c('spuzzum_creek','upper_squamish','chilliwack_river','tenquille_lake')
site.names <- c('Spuzzum Creek','Upper Squamish','Chilliwack River','Tenquille Lake')
slen <- 100

swe.sims <- matrix(0,nrow=slen,ncol=13696)
snow.sims <- matrix(0,nrow=slen,ncol=13696)


    plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
    png(file=paste0(plot.dir,model,'.snow.pillow.series.png'),width=1000,height=900)
    par(mfrow=c(4,1),oma=c(1,2,1,1))    


##Loop over sites

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)
    ##Reanalysis 800m data
    if (site=='spuzzum_creek') {
        clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sites/sp_testing/',site,'_',model,'_800m_data_193_121.csv')
    } else {
        clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sites/',site,'_',model,'_800m_data.csv')
    }
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

    coeffs <- list(a=-49.49,b=0.4128,c=2.6545,d=1.0209)
    for (k in 1:slen) {
    print(k)
    test.snow <- hyper.snow(pr.data,tasmax.data,tasmin.data,coeffs)
    results <- snow.melt(Date=dates, precip_mm=pr.data, Tmax_C=tasmax.data, Tmin_C=tasmin.data,Snow=test.snow,
                         lat_deg=lat.bnds, slope=site.slope, aspect=site.aspect, tempHt=1, windHt=2, groundAlbedo=0.25,
                         SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)
    swe.sims[k,] <- results$swe
    snow.sims[k,] <- results$snowdepth
    }                                      

    winter.flag <- grep('(*-12-*|*-01-*|*-02-*)',dates)
    swe.winter <- swe.sims[,winter.flag]*1000
    swe.winter.pos <- swe.winter[swe.winter>0]

    pack.winter <- snow.sims[,winter.flag]*100
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
    sb <- 1:length(pillow.dates)
    date.subset <- format(as.Date(dates),'%Y-%m-%d') %in% pillow.dates[sb]
    pillow.subset <- pillow.dates %in% format(as.Date(dates),'%Y-%m-%d')

 
    pillow.winter.flag <- grep('(*-12-*|*-01-*|*-02-*)',pillow.dates)
    pillow.swe.winter <- pillow.swe[pillow.winter.flag]    
    p.swe.winter.pos <- pillow.swe.winter[pillow.swe.winter>0]

    pillow.pack.winter <- pillow.pack[pillow.winter.flag]    
    p.pack.winter.pos <- pillow.pack.winter[pillow.pack.winter>0]


    ##SWE
    if (type=='SWE') {
      par(mar=c(4,6,2,2))
      plot(as.Date(dates[date.subset]),round(results$swe[date.subset]*1000,2),xlab='Date',ylab='SWE (mm)',yaxs='i',
           type='l',lwd=3,col='blue',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
           xlim=c(as.Date('1992-06-01'),as.Date('2011-01-01')),ylim=c(0,3000))
      apply(swe.sims[,date.subset],1,function(x,y){lines(y,x*1000,col='lightblue',lwd=2.5)},as.Date(dates[date.subset]))
      points(as.Date(pillow.dates),pillow.swe,pch=16,col='black')
      lines(as.Date(dates[date.subset]),apply(swe.sims*1000,2,mean)[date.subset],col='blue',lwd=3.5)
      text(as.Date('1994-01-01'),2500,site.names[i],cex=2.5)
      if (i==4) { 
         legend('bottomleft',legend=c('ASP Obs.','Model','Model Mean'),col=c('black','lightblue','blue'),pch=16,cex=1.75)
      }
    }

    ##Snowpack                   
    if (type=='pack') {                
      plot(as.Date(dates)[date.subset],results$snowdepth[date.subset]*100,type='l',lwd=3,col='blue',cex.axis=1.5,yaxs='i',
               xlab='Date',ylab='Snowpack (cm)',main=site.names[i], cex.lab=1.75,cex.axis=1.75,cex.main=2)
      apply(snow.sims[,date.subset],1,function(x,y){lines(y,x*100,col='lightblue',lwd=2.5)},as.Date(dates[date.subset]))
      points(as.Date(pillow.dates),pillow.pack,cex=1,col='black',pch=16)
      lines(as.Date(dates[date.subset]),apply(snow.sims*100,2,mean)[date.subset],col='blue',lwd=3.5)
      abline(h=0)
    }
}


    dev.off()