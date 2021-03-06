##SCript to compute the start and end dates and lengths of the snow season for each year

library(ncdf4)
library(plotrix)

get.coordinates <- function(site) {

  coordinates <- list(spuzzum_creek=c(-121.686,49.674,1197),
                      chilliwack_river=c(-121.71667,49.0333,1600),
                      upper_squamish=c(-123.4333,50.1500,1340),
                      tenquille_lake=c(-122.9333,50.5333,1680))
  rv <- coordinates[[site]]
  return(rv)
}

compare.precip <- function(sim.dates,pillow.dates,season.dates,sim.pr,pillow.pr) {
   
   clim.match <- rep(0,nrow(season.dates))
   pillow.match <- rep(0,nrow(season.dates))

   for (j in 1:nrow(season.dates)) {
     sim.bnds <- c(grep(season.dates$Start[j],sim.dates),grep(season.dates$End[j],sim.dates))
     asp.bnds <- c(grep(season.dates$Start[j],pillow.dates),grep(season.dates$End[j],pillow.dates))
     clim.match[j] <- sum(sim.pr[sim.bnds[1]:sim.bnds[2]],na.rm=T)
     pillow.match[j] <- sum(pillow.pr[asp.bnds[1]:asp.bnds[2]],na.rm=T)
     print(sum(is.na(pillow.pr[asp.bnds[1]:asp.bnds[2]])))
     print(diff(asp.bnds))
   }
}


##Slope and Aspect Values
model <- 'ERA'
type <- 'SWE'


sites <- c('spuzzum_creek','upper_squamish','chilliwack_river','tenquille_lake')
site.names <- c('Spuzzum Creek','Upper Squamish','Chilliwack River','Tenquille Lake')

slen <- 101

plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
png(file=paste0(plot.dir,'ERA.PNWNAmet.SWE.pillow.accumulated.precipitation.png'),width=1000,height=900)
par(mfrow=c(4,1),oma=c(1,2,1,1))    

results.matrix <- matrix(0,nrow=15,ncol=length(sites))

##Loop over sites

for (i in seq_along(sites)) {
    
    site <- sites[i]
    print(site)
    clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_',model,'_800m_data.csv')
    clim.data <- read.csv(clim.file,header=T,as.is=T)
    sim.dates <- clim.data$Dates

    pnw.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_PNWNAmet_800m_data.csv')
    pnw.data <- read.csv(pnw.file,header=T,as.is=T)
    pnw.dates <- pnw.data$Dates

    dates.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/snow_seasons/',site,'_season_dates.csv') 
    season.dates <- read.csv(dates.file,header=T,as.is=T)
    season.dates <- season.dates[as.Date(season.dates$End) < as.Date('2012-12-31'),]
    pillow.starts <- as.numeric(format(as.Date(season.dates$Start),'%j'))
    pillow.ends <- as.numeric(format(as.Date(season.dates$End),'%j'))
    pillow.lengths <- as.numeric(season.dates$Length)

    coords <- get.coordinates(site)
    lat.bnds <- coords[2]
    elev <- coords[3]

    ##Snow Pillow Data
    ##pillow.file <- paste('/storage/data/projects/rci/data/assessments/snow_model/snow_pillow/',site,'_asp.csv',sep='')
    pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'_asp.csv',sep='')
    pillow.data <- read.csv(pillow.file,header=T,as.is=T)
    pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
    pillow.tasmax <- pillow.data[,3]
    pillow.tasmin <- pillow.data[,5]
    pillow.tas <- (pillow.tasmax + pillow.tasmin)/2
    pillow.precip <- pillow.data[,7]##mm
    pillow.swe <- pillow.data[,11] ##mm
    pillow.pack <- pillow.data[,13] ##cm

    date.subset <- format(as.Date(sim.dates),'%Y-%m-%d') %in% pillow.dates
    pillow.subset <- pillow.dates %in% format(as.Date(sim.dates),'%Y-%m-%d')
    
    cumulative.sim.pr <- clim.data$Pr*0
    cumulative.pnw.pr <- pnw.data$Pr*0
    cumulative.asp.pr <- pillow.precip*0

    se.len <- nrow(season.dates)
    
    for (j in 1:se.len) {

      sst <- grep(season.dates$Start[j],sim.dates)  
      sen <- grep(season.dates$Peak[j],sim.dates)  
      wst <- grep(season.dates$Start[j],pnw.dates)  
      wen <- grep(season.dates$Peak[j],pnw.dates)  
      pst <- grep(season.dates$Start[j],pillow.dates)  
      pen <- grep(season.dates$End[j],pillow.dates)  
      pslct <- (pillow.dates < pillow.dates[pen]) & (pillow.dates > pillow.dates[pst])       
      cumulative.sim.pr[sst:sen] <- cumsum(clim.data$Pr[sst:sen])    
      cumulative.pnw.pr[wst:wen] <- cumsum(pnw.data$Pr[wst:wen])    
      cumulative.asp.pr[pst:pen] <- pillow.swe[pst:pen] ##   cumsum(pillow.precip[pst:pen])    


    }
    ##plot(as.Date(sim.dates),cumulative.sim.pr,xlim=as.Date(c('1990-01-01','2012-12-31')))
    par(mar=c(4,6,2,2))

    plot(as.Date(sim.dates),cumulative.sim.pr,xlim=as.Date(c('1992-06-01','2012-08-31')),
    col='blue',type='l',lwd=3,ylim=c(0,3000),
    xlab='Date',ylab='SWE (mm)',yaxs='i',
    cex.axis=1.75,cex.lab=1.75,cex.main=2)
    text(as.Date('1996-01-01'),2500,site.names[i],cex=2.5)
    lines(as.Date(pnw.dates),cumulative.pnw.pr,xlim=as.Date(c('1990-01-01','2012-12-31')),col='green',lwd=3)
    points(as.Date(pillow.dates),cumulative.asp.pr,pch=18,col='black',cex=1.25)
    
##    plot(as.Date(sim.dates)[date.subset],clim.data$Tas[date.subset]-pillow.tas[pillow.subset],xlim=as.Date(c('1998-01-01','2018-12-31')),col='blue',type='l')
    ##lines(as.Date(pillow.dates),pillow.tas,lwd=3,col='black')
    print(mean(clim.data$Tas[date.subset]-pillow.tas[pillow.subset],na.rm=T))
    if (i==4) {
       legend('bottomleft',legend=c('ASP Obs.','ERA-Int','PNWNAmet'),col=c('black','blue','green'),pch=16,cex=1.75)
    }
          


}
dev.off()

browser()

      par(mar=c(4,6,2,2))
      plot(as.Date(sim.dates[date.subset]),apply(swe.sims[date.subset,],1,mean),xlab='Date',ylab='SWE (mm)',yaxs='i',
           type='l',lwd=3,col='blue',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
           xlim=c(as.Date('1998-06-01'),as.Date('2002-01-01')),ylim=c(0,3000))
      apply(swe.sims[date.subset,],2,function(x,y){lines(y,x,col='lightblue',lwd=2.5)},as.Date(sim.dates[date.subset]))
      points(as.Date(pillow.dates),pillow.swe,pch=16,col='black')
      lines(as.Date(sim.dates[date.subset]),apply(swe.sims,1,mean)[date.subset],col='blue',lwd=3.5)
##      text(as.Date('1994-01-01'),2500,site.names[i],cex=2.5)
##      if (i==4) { 
##         legend('bottomleft',legend=c('ASP Obs.','Model','Model Mean'),col=c('black','lightblue','blue'),pch=16,cex=1.75)
##      }
    

      par(mar=c(4,6,2,2))
      plot(as.Date(sim.dates[date.subset]),clim.data$Pr,xlab='Date',ylab='SWE (mm)',yaxs='i',
           type='l',lwd=3,col='blue',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
           xlim=c(as.Date('1998-06-01'),as.Date('2002-01-01')),ylim=c(0,3000))
      apply(swe.sims[date.subset,],2,function(x,y){lines(y,x,col='lightblue',lwd=2.5)},as.Date(sim.dates[date.subset]))
      points(as.Date(pillow.dates),pillow.swe,pch=16,col='black')
      lines(as.Date(sim.dates[date.subset]),apply(swe.sims,1,mean)[date.subset],col='blue',lwd=3.5)
      



if (1==0) {
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


##    dev.off()