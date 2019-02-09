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

##Slope and Aspect Values
model <- 'ERA'
type <- 'SWE'


sites <- c('spuzzum_creek','upper_squamish','chilliwack_river','tenquille_lake')
site.names <- c('Spuzzum Creek','Upper Squamish','Chilliwack River','Tenquille Lake')
slen <- 100

    plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
    png(file=paste0(plot.dir,model,'.SWE.normalized.pillow.series.1001.2018.png'),width=1000,height=900)
    par(mfrow=c(4,1),oma=c(1,2,1,1))    


##Loop over sites

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)
    ##if (site=='spuzzum_creek') {
    ##    clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sites/sp_testing/',site,'_',model,'_800m_data_193_121.csv')
    ##} else {
        clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_',model,'_800m_data.csv')
    ##}
    clim.data <- read.csv(clim.file,header=T,as.is=T)
    dates <- clim.data$Dates

    pack.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sims/',site,'.',tolower(model),'.snow.1001.csv')
    pack.sims <- read.csv(pack.file,header=T,as.is=T)
    swe.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sims/',site,'.',tolower(model),'.swe.1001.csv')
    swe.sims <- as.matrix(read.csv(swe.file,header=T,as.is=T))

    coords <- get.coordinates(site)
    lat.bnds <- coords[2]
    elev <- coords[3]

    winter.flag <- grep('(*-12-*|*-01-*|*-02-*)',dates)
    swe.winter <- swe.sims[winter.flag,]
    swe.winter.pos <- swe.winter[swe.winter>0]

    pack.winter <- pack.sims[winter.flag,]
    pack.winter.pos <- pack.winter[pack.winter>0]

    ##Snow Pillow Data
    pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'_asp.csv',sep='')
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


    sims <- apply(swe.sims,1,mean)[date.subset]
    asp <- pillow.swe[pillow.subset]
    flag <- is.na(asp)

    sims.na <- sims[!flag]
    asp.na <- asp[!flag]
    print(cor(sims.na,asp.na))


    ##Normalized Series
    if (1==0) {
      par(mar=c(4,6,2,2))
      swe.sim.means <- apply(swe.sims[date.subset,],2,mean)
      swe.sim.sds <- apply(swe.sims[date.subset,],2,sd)
      swe.sim.norms <- (swe.sims - swe.sim.means)/swe.sim.sds      
      
      asp.mean <- mean(pillow.swe,na.rm=T)
      asp.sd   <- sd(pillow.swe,na.rm=T)
      pillow.swe.norm <- (pillow.swe - asp.mean)/asp.sd

      plot(as.Date(dates[date.subset]),apply(swe.sim.norms[date.subset,],1,mean),xlab='Date',ylab='SWE (mm)',yaxs='i',
           type='l',lwd=3,col='blue',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
           xlim=c(as.Date('1992-06-01'),as.Date('2019-01-01')),ylim=c(-2,6))
      apply(swe.sim.norms[date.subset,],2,function(x,y){lines(y,x,col='lightblue',lwd=2.5)},as.Date(dates[date.subset]))
      points(as.Date(pillow.dates),pillow.swe.norm,pch=16,col='black')
      lines(as.Date(dates[date.subset]),apply(swe.sim.norms,1,mean)[date.subset],col='blue',lwd=3.5)
      ##text(as.Date('1994-01-01'),2500,site.names[i],cex=2.5)
      if (i==4) { 
         legend('bottomleft',legend=c('ASP Obs.','Model','Model Mean'),col=c('black','lightblue','blue'),pch=16,cex=1.75)
      }
    }
      
    if (1==1) {
    ##SWE
    if (type=='SWE') {
      par(mar=c(4,6,2,2))
      plot(as.Date(dates[date.subset]),apply(swe.sims[date.subset,],1,mean),xlab='Date',ylab='SWE (mm)',yaxs='i',
           type='l',lwd=3,col='blue',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,xaxs='i',
           xlim=c(as.Date('1992-08-01'),as.Date('2018-08-01')),ylim=c(0,3000))
      apply(swe.sims[date.subset,],2,function(x,y){lines(y,x,col='lightblue',lwd=2.5)},as.Date(dates[date.subset]))
      points(as.Date(pillow.dates),pillow.swe,pch=16,col='black')
      lines(as.Date(dates[date.subset]),apply(swe.sims,1,mean)[date.subset],col='blue',lwd=3.5)
      text(as.Date('1995-01-01'),2500,site.names[i],cex=2.5)
      if (i==4) { 
         legend('bottomleft',legend=c('ASP Obs.',model,paste0(model,' Mean')),col=c('black','lightblue','blue'),pch=16,cex=1.75)
      }
    }

    ##Snowpack                   
    if (type=='pack') {                
      plot(as.Date(dates)[date.subset],apply(pack.sims[date.subset,],1,mean),type='l',lwd=3,col='white',cex.axis=1.5,yaxs='i',
               xlab='Date',ylab='Snowpack (cm)',main=site.names[i], cex.lab=1.75,cex.axis=1.75,cex.main=2)
      apply(pack.sims[date.subset,],2,function(x,y){lines(y,x,col='lightblue',lwd=2.5)},as.Date(dates[date.subset]))
      points(as.Date(pillow.dates),pillow.pack,cex=1,col='black',pch=16)
      lines(as.Date(dates[date.subset]),apply(pack.sims,1,mean)[date.subset],col='blue',lwd=3.5)
      abline(h=0)
    }
    }
}


    dev.off()