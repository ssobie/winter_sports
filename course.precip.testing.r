##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)


hyper_snow <- function(pr,tasmax,tasmin) {
           
    coeffs <- list(a=-49.49,b=0.5628,c=1.5,d=1.0209)

    tas <- (tasmax+tasmin)/2                # degrees C
    frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)
    sample <- runif(length(tas),min=0,max=100)
    test <- sample > frac
    snow.type <- rep(TRUE,length(tas))
    snow.type[test] <- FALSE

    NewSnowWatEq <- pr
    NewSnowWatEq[!snow.type] <- 0
    R_m <- pr
    R_m[snow.type] <- 0
    rv <- list(snow=NewSnowWatEq,
               rain=R_m)
    return(NewSnowWatEq)
}



get_snow_season_dates <- function(input.snow,dates) {

  years <- format(as.Date(dates),'%Y')
  jdays <- format(as.Date(dates),'%j')

  peaks.ix <- tapply(input.snow,as.factor(years),function(x){which.max(x[1:300])})
  peaks <- as.Date(paste0(levels(as.factor(years)),'-',sprintf('%03d',peaks.ix)),format='%Y-%j')
  
##Fix to return the date of peak SWE
  sim.snow <- input.snow
  sim.snow[sim.snow >= 25] <- 100
  sim.snow[sim.snow < 25] <- 0

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
  len <- length(ix)
  lengths <- ix[seq(2,len,2)] - ix[seq(1,len-1,2)]
  ###starts <- as.numeric(format(as.Date(dates[ix[seq(1,len-1,2)]]),'%j'))
  ###ends <- as.numeric(format(as.Date(dates[ix[seq(2,len,2)]]),'%j'))
  starts <- as.Date(dates[ix[seq(1,len-1,2)]])
  ends <- as.Date(dates[ix[seq(2,len,2)]])

  peaks.match <- as.Date('1900-01-01')
  for (k in seq_along(starts)) {
     peaks.match <- c(peaks.match,peaks[which.min(abs(peaks-starts[k]))])
  }

  rv <- list(lengths=lengths,
             starts =starts,
             peaks = peaks.match[-1],
             ends = ends,
             years = unique(years))
  return(rv)
}




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


sites <- c('callaghan',
           'palisade_lake',
           'stave_lake',
           'nahatlatch',
           'klesilkwa',
           'brookmere')

site.names <- c('Callaghan','Palisade Lake',
                'Stave Lake','Nahatlatch',
                'Klesilkwa','Brookmere')


sites <- c('orchid_lake',           
           'grouse_mountain',
           'dog_mountain',                      
           'wahleach',           
           'lightning_lake',           
           'shovelnose_mountain',
           'hamilton_hill')

site.names <- c('Orchid Lake','Grouse Mountain',
                'Dog Mountain','Wahleach',
                'Lightning Lake','Shovelnose Mountain',
                'Hamilton Hill')



##North Eastern Sites
sites <- c('brookmere',
           'shovelnose_mountain',
           'hamilton_hill')
           
site.names <- c('Brookmere',
                'Shovelnose Mountain',
                'Hamilton Hill')


##Selected
sites <- c('grouse_mountain',
           'nahatlatch',
           'brookmere')

site.names <- c('Grouse Mountain',
           'Nahatlatch',
           'Brookmere')

##Central Sites
sites <- c('callaghan',
           'stave_lake',
           'nahatlatch')

site.names <- c('Callaghan',
                'Stave Lake',
                'Nahatlatch')

##South Eastern Sites                
sites <- c('lightning_lake',
           'klesilkwa',
           'wahleach')

site.names <- c('Lightning Lake',
                'Klesilkwa',
                'Wahleach')

##Eastern Sites
sites <- c('lightning_lake',
           'brookmere',
           'shovelnose_mountain',
           'hamilton_hill')
site.names <- c('Lightning Lake',
                'Brookmere',
                'Shovelnose Mountain',
                'Hamilton Hill')

##North Shore Sites
sites <- c('grouse_mountain',
           'dog_mountain',
           'orchid_lake',           
           'palisade_lake')                      
site.names <- c('Grouse Mountain','Dog Mountain',
                'Orchid Lake','Palisade Lake')



model <- 'ERA'

course.site.swe <- vector(mode='list',length=length(sites))
model.site.swe <- vector(mode='list',length=length(sites))


model.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/grouse_mountain_',model,'_800m_data.csv')
model.data <- read.csv(model.file,header=T,as.is=T)


##Loop over sites
model.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sims/'
plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'

#png(file=paste0(plot.dir,model,'course.precip.compare.eastern.sites.png'),width=7,height=6,units='in',res=600,pointsize=6,bg='white')
#par(mfrow=c(4,1))    

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)
    ##Reanalysis 800m data
    clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_',model,'_800m_data.csv')
    clim.data <- read.csv(clim.file,header=T,as.is=T)
    clim.pr <- clim.data$Pr
    clim.tasmax <- clim.data$Tasmax
    clim.tasmin <- clim.data$Tasmin
    tas.data <- clim.data$Tas
    clim.dates <- clim.data$Dates

    ##Apply snow fall model
    slen <- 10
    snow.series <- matrix(0,nrow=length(clim.pr),ncol=slen)
    for (k in 1:slen) {
       snow.series[,k] <- hyper_snow(clim.pr,clim.tasmax,clim.tasmin)
    }

    pnw.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_PNWNAmet_800m_data.csv')
    pnw.data <- read.csv(pnw.file,header=T,as.is=T)
    pnw.dates <- pnw.data$Dates
    pnw.pr <- pnw.data$Pr

    coords <- get.coordinates(site)
    lat.bnds <- coords[2]
    elev <- coords[3]

    ##Snow Course Data
    course.file <- paste0('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
    ##course.file <- paste('/storage/data/projects/rci/data/assessments/snow_model)
    course.data <- read.csv(course.file,header=T,as.is=T)
    course.dates <- format(as.Date(course.data[,1]),'%Y-%m-%d')
    course.swe <- course.data[,3] ##mm

    swe.sims <- read.csv(paste0(model.dir,site,'.',tolower(model),'.swe.1001.csv'),header=T,as.is=T)
    swe.mean <- apply(swe.sims,1,mean,na.rm=T)
    rm(swe.sims)

    season.dates <- get_snow_season_dates(swe.mean,clim.dates)
    se.len <- length(season.dates$lengths)
    print(season.dates$starts)
    cumulative.sim.pr <- clim.data$Pr*0
    cumulative.sim.swe <- snow.series*0
    cumulative.pnw.pr <- pnw.data$Pr*0

    for (j in 1:se.len) {

      sst <- grep(season.dates$starts[j],clim.dates)
      sen <- grep(season.dates$peaks[j],clim.dates)

##      wst <- grep(as.Date(paste0(season.dates$years[j-1],'-',sprintf('%03d',season.dates$starts[j])),format='%Y-%j'),pnw.dates)
##      wen <- grep(as.Date(paste0(season.dates$years[j],'-',sprintf('%03d',season.dates$peaks[j])),format='%Y-%j'),pnw.dates)
      cumulative.sim.pr[sst:sen] <- cumsum(clim.data$Pr[sst:sen])
      for (k in 1:slen) {
        cumulative.sim.swe[sst:sen,k] <- cumsum(snow.series[sst:sen,k])
      }
##      cumulative.pnw.pr[wst:wen] <- cumsum(pnw.data$Pr[wst:wen])
    }

##Snow Course Comparison
    yupp <- max(c(max(course.swe,na.rm=T),max(cumulative.sim.pr),max(cumulative.pnw.pr)))
    ymax <- max(c(max(course.swe,na.rm=T),max(cumulative.sim.pr),max(cumulative.pnw.pr)))

    par(mar=c(5.1,5,2.1,2.1))
    plot(as.Date(course.dates),course.swe,cex=1.5,col='black',pch=16,
             xlim=c(as.Date('1999-09-01'),as.Date('2000-04-30')),ylim=c(0,ymax),
             main='',xlab='Date',ylab='SWE (mm)', cex.lab=1.95,axes=F)
    axis(1,at=as.Date(c('1998-01-01','2000-01-01','2002-01-01','2004-01-01')),label=c('1998','2000','2002','2004'),cex.axis=1.95)
    axis(2,at=seq(0,yupp,round((yupp-100)/3,-2)),label=seq(0,yupp,round((yupp-100)/3,-2)),cex.axis=1.95)

    ##lines(as.Date(clim.dates),swe.mean,lwd=3,col='blue')
    apply(cumulative.sim.swe,2,function(y,x){lines(x,y,col='gray',lwd=1)},as.Date(clim.dates))
    lines(as.Date(clim.dates),apply(cumulative.sim.swe,1,mean),col='darkgreen',lwd=2)
    lines(as.Date(clim.dates),cumulative.sim.pr,xlim=as.Date(c('1990-01-01','2012-12-31')),col='blue',lwd=2)

##    lines(as.Date(pnw.dates),cumulative.pnw.pr,xlim=as.Date(c('1990-01-01','2012-12-31')),col='green',lwd=3)

    points(as.Date(course.dates),course.swe,cex=1.75,col='black',pch=16)
    text(as.Date('2005-01-01'),0.95*ymax,site.names[i],cex=2)
    if (i==4) {
       legend('topleft',legend=c('Course Obs.','ERA','PNWNAmet'),col=c('black','blue','green'),pch=16,cex=1.75)
    }
    box(which='plot')
    abline(h=0)

browser()
}        

##dev.off()

