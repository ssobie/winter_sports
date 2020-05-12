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

get.snow.season.dates <- function(input.snow,dates) {

  years <- format(as.Date(dates),'%Y')
  jdays <- format(as.Date(dates),'%j')
  jdays.list <- tapply(jdays,as.factor(years),list)
  peaks <- tapply(input.snow,as.factor(years),function(x){which.max(x[1:300])})                      
                    
  sim.snow <- input.snow
  sim.snow[sim.snow>0] <- 100
                      
  patterns <- rle(sim.snow)
  pt     <- patterns$lengths
  pt.sum <- cumsum(patterns$lengths)
  all.zeroes <- which(patterns$values==100)
  zeroes <- all.zeroes[patterns$lengths[all.zeroes] > 100]
  top.ix <- pt.sum[zeroes]
  bot.ix <- pt[zeroes]

  pt.diff <- c(pt[1],top.ix[1], top.ix[-1] - bot.ix[-1], tail(top.ix,1))
  ix <- sort(unique(c(pt.diff,top.ix))) 
  if (ix[1] > 100) {
     ix <- sort(unique(c(1,pt.diff,top.ix)))
  }   
##  plot(as.Date(dates),input.snow)
##  abline(v=as.Date(dates[ix]),col='red')
  len <- length(ix)

  lengths <- ix[seq(2,len,2)] - ix[seq(1,len-1,2)]
  starts <- as.numeric(format(as.Date(dates[ix[seq(1,len-1,2)]]),'%j'))
  ends <- as.numeric(format(as.Date(dates[ix[seq(2,len,2)]]),'%j'))

  rv <- list(lengths=lengths,
             starts =starts,
             peaks = peaks,
             ends = ends)

  return(rv)               

}

##Slope and Aspect Values
model <- 'NCEP2'
type <- 'SWE'


sites <- c('spuzzum_creek','upper_squamish','chilliwack_river','tenquille_lake')
site.names <- c('Spuzzum Creek','Upper Squamish','Chilliwack River','Tenquille Lake')
slen <- 1001

results.matrix <- matrix(0,nrow=20,ncol=length(sites))

##Loop over sites

for (i in seq_along(sites)) {

    site <- sites[i]
    print(site)
    clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_',model,'_800m_data.csv')

    clim.data <- read.csv(clim.file,header=T,as.is=T)
    sim.dates <- clim.data$Dates

    dates.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/snow_seasons/',site,'_season_dates.csv') 
    season.dates <- read.csv(dates.file,header=T,as.is=T)
    pillow.starts <- as.numeric(format(as.Date(season.dates$Start),'%j'))
    pillow.peaks <- as.numeric(format(as.Date(season.dates$Peak),'%j'))
    pillow.ends <- as.numeric(format(as.Date(season.dates$End),'%j'))
    pillow.lengths <- as.numeric(season.dates$Length)

    pack.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sims/',site,'.',tolower(model),'.snow.1001.csv')
    pack.sims <- read.csv(pack.file,header=T,as.is=T)
    swe.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sims/',site,'.',tolower(model),'.swe.1001.csv')
    swe.sims <- read.csv(swe.file,header=T,as.is=T)

    sim.years <- unique(format(as.Date(sim.dates),'%Y'))
    ylen <-  length(sim.years)
    test.dates <- get.snow.season.dates(swe.sims[,1],sim.dates)
    sim.lengths <- matrix(0,nrow=slen,ncol=ylen)
    sim.starts <- matrix(0,nrow=slen,ncol=ylen)
    sim.peaks <- matrix(0,nrow=slen,ncol=ylen)
    sim.ends <- matrix(0,nrow=slen,ncol=ylen)
    
    for (s in 1:slen) {
      swe.se.dates <- get.snow.season.dates(swe.sims[,s],sim.dates)
      sim.lengths[s,1:length(swe.se.dates$lengths)] <- swe.se.dates$lengths
      sim.starts[s,1:length(swe.se.dates$starts)] <- swe.se.dates$starts
      sim.peaks[s,1:length(swe.se.dates$peaks)] <- swe.se.dates$peaks
      sim.ends[s,1:length(swe.se.dates$ends)] <- swe.se.dates$ends
    }
    sim.lengths[sim.lengths==0] <- NA
    sim.starts[sim.starts==0] <- NA
    sim.peaks[sim.peaks==0] <- NA
    sim.ends[sim.ends==0] <- NA

    year.match <- as.numeric(sim.years) %in% as.numeric((format(as.Date(season.dates$Start),'%Y')))
    
    lengths.diff <- season.dates$Length - sim.lengths[,year.match]
    starts.diff <-  pillow.starts - sim.starts[,year.match]
    peaks.diff <-  pillow.peaks - sim.peaks[,year.match]
    ends.diff <-  pillow.ends - sim.ends[,year.match]

    rv <- c(mean(pillow.starts),mean(pillow.peaks),mean(pillow.ends),mean(pillow.lengths),
            mean(sim.starts[,year.match]),mean(sim.peaks[,year.match]),mean(sim.ends[,year.match]),mean(sim.lengths[,year.match]),
            sd(sim.starts[,year.match]),sd(sim.peaks[,year.match]),sd(sim.ends[,year.match]),sd(sim.lengths[,year.match]),
            mean(starts.diff),mean(peaks.diff),mean(ends.diff),mean(lengths.diff),
            sd(starts.diff),sd(peaks.diff),sd(ends.diff),sd(lengths.diff))
    results.matrix[,i] <- round(rv,1)

}

    results.matrix <- cbind(c('Site','ASP Start','ASP Peak','ASP End','ASP Length',
                              'Sim Start','Sim Peak','Sim End','Sim Length',
                              'Start SD' ,'Peak SD','End SD' ,'Length SD' ,
                              'Start Diff','Peak Diff','End Diff','Length Diff',      
                              'SDiff SD'  ,'PDiff SD','EDiff SD','LDiff SD'),rbind(site.names,results.matrix))



