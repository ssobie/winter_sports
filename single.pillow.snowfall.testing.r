##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/snow.model.functional.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)

    ##coeffs <- list(a=-51.8,b=0.875,c=4.1,d=1.0209) ##Blackwall Pea
    ##coeffs <- list(a=-49.45,b=0.5,c=2.25,d=1.0209) ##Tenquille Lake
    ##coeffs <- list(a=-51.8,b=0.875,c=4.1,d=1.0209) ##Spuzzum Cree
    ##coeffs <- list(a=-51.8,b=0.9,c=4.0,d=1.0209) ##Chilliwack River
    ##coeffs <- list(a=-51.7,b=0.575,c=3.7,d=1.0209) ##Upper Squamish

    coeffs <- list(a=-51.45,b=1.025,c=3.975,d=1.0209) ##Upper Squamish
    coeffs <- list(a=-51.9,b=0.95,c=4.5,d=1.0209) ##Chilliwack River
    coeffs <- list(a=-51.7,b=0.625,c=4.5,d=1.0209) ##Spuzzum Creek
    coeffs <- list(a=-48.65,b=0.75,c=1.575,d=1.0209) ##Tenquille Lake
    coeffs <- list(a=-48.0,b=0.25,c=2.675,d=1.0209) ##Blackwall Peak
    coeffs <- list(a=-51.85,b=0.8,c=4.25,d=1.0209) ##Wahleach Lake





hyper_snow <- function(pr,tasmax,tasmin) {
           
    ##coeffs <- list(a=-48.49,b=0.7628,c=3.0,d=1.0209) ##Best manual guess



    tas <- (tasmax+tasmin)/2                # degrees C
    frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)
    sample <- runif(length(tas),min=0,max=100)
    test <- sample > frac
    high.temp <- tas > 15
    test[high.temp] <- TRUE
    low.temp <- tas < -10
    test[low.temp] <- FALSE

    snow.type <- rep(TRUE,length(tas))
    snow.type[test] <- FALSE

    NewSnowWatEq <- pr
    NewSnowWatEq[!snow.type] <- 0
    R_m <- pr
    R_m[snow.type] <- 0
    rv <- list(swe=NewSnowWatEq,
               rain=R_m)
    return(rv)
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

##North Shore Sites
site <- 'wahleach_lake'
site.name <- 'Wahleach'

model <- 'PNWNAmet'

##Loop over sites
model.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sims/'

print(site)
##Reanalysis 800m data
clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_',model,'_800m_data.csv')
clim.data <- read.csv(clim.file,header=T,as.is=T)
clim.pr <- clim.data$Pr
clim.tasmax <- clim.data$Tasmax
clim.tasmin <- clim.data$Tasmin
clim.tas <- clim.data$Tas
clim.dates <- clim.data$Dates

##Apply snow fall model
slen <- 10
snow.series <- matrix(0,nrow=length(clim.pr),ncol=slen)
rain.series <- matrix(0,nrow=length(clim.pr),ncol=slen)
for (k in 1:slen) {
   snow <- hyper_snow(clim.pr,clim.tasmax,clim.tasmin)
   snow.series[,k] <- snow$swe
   rain.series[,k] <- snow$rain
}

snow.series.mean <- apply(snow.series,1,mean)
rain.series.mean <- apply(snow.series,1,mean)
rain.series.mean[snow.series.mean !=0] <- 0 

coords <- get_coordinates(site)
lat.bnds <- coords[2]
elev <- coords[3]

##Snow Pillow Data
pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'.csv',sep='')
pillow.data <- read.csv(pillow.file,header=T,as.is=T)
pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
pillow.swe <- pillow.data[,11] ##mm

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

site.slope <- bc.slopes[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]
site.aspect <- bc.aspects[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]

snow.mean <- list(swe=snow.series.mean*0.001,rain=rain.series.mean*0.001)
results <- snow_melt_functional(Date=clim.dates, precip_mm=snow.mean, Tmax_C=clim.tasmax, Tmin_C=clim.tasmin,
                     lat_deg=coords[2], slope=site.slope, aspect=site.aspect, tempHt=1, windHt=2, groundAlbedo=0.25,
                      SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)

season.dates <- get_snow_season_dates(results$swe*1000,clim.dates)
se.len <- length(season.dates$lengths)
print(season.dates$starts)
cumulative.sim.pr <- clim.data$Pr*0
cumulative.end.pr <- clim.data$Pr*0
cumulative.snow.mean <- clim.data$Pr*0
cumulative.snow.end <- clim.data$Pr*0

cumulative.sim.swe <- snow.series*0
cumulative.end.swe <- snow.series*0

for (j in 1:se.len) {
   sst <- grep(season.dates$starts[j],clim.dates)
   sen <- grep(season.dates$peaks[j],clim.dates)
   een <- grep(season.dates$ends[j],clim.dates)

   cumulative.sim.pr[sst:sen] <- cumsum(clim.data$Pr[sst:sen])
   cumulative.end.pr[sst:een] <- cumsum(clim.data$Pr[sst:een])

   cumulative.snow.mean[sst:sen] <- cumsum(snow.series.mean[sst:sen])
   cumulative.snow.end[sst:een] <- cumsum(snow.series.mean[sst:een])

   for (k in 1:slen) {
     cumulative.sim.swe[sst:sen,k] <- cumsum(snow.series[sst:sen,k])
     cumulative.end.swe[sst:een,k] <- cumsum(snow.series[sst:een,k])
   }
}

##Snow Course Comparison
yupp <- max(c(max(pillow.swe,na.rm=T),max(cumulative.sim.pr)))
ymax <- max(c(max(pillow.swe,na.rm=T),max(cumulative.sim.pr)))


dst <- '2005-10-01'
den <- '2008-06-30'

par(mfrow=c(2,1))
par(mar=c(4.1,5,1.1,1.1))
if (1==0) {
plot(as.Date(pillow.dates),pillow.swe,cex=1.5,col='black',pch=16,
     xlim=c(as.Date(dst),as.Date(den)),ylim=c(0,ymax),
     main='',xlab='Date',ylab='SWE (mm)', cex.lab=1.25,axes=F)
axis(1,at=as.Date(c('1998-01-01','2000-01-01','2002-01-01','2004-01-01')),label=c('1998','2000','2002','2004'),cex.axis=1.25)
axis(2,at=seq(0,yupp,round((yupp-100)/3,-2)),label=seq(0,yupp,round((yupp-100)/3,-2)),cex.axis=1.25)

##lines(as.Date(clim.dates),swe.mean,lwd=3,col='blue')
##apply(cumulative.end.swe,2,function(y,x){lines(x,y,col='gray',lwd=1)},as.Date(clim.dates))
##lines(as.Date(clim.dates),apply(cumulative.end.swe,1,mean),col='lightblue',lwd=2)

##apply(cumulative.sim.swe,2,function(y,x){lines(x,y,col='gray',lwd=1)},as.Date(clim.dates))
##lines(as.Date(clim.dates),apply(cumulative.sim.swe,1,mean),col='darkgreen',lwd=2)
##lines(as.Date(clim.dates),cumulative.sim.pr,col='blue',lwd=2)
##lines(as.Date(clim.dates),cumulative.end.pr,col='blue',lwd=2,lty=2)

lines(as.Date(clim.dates),cumulative.snow.mean,col='blue',lwd=2)
lines(as.Date(clim.dates),cumulative.snow.end,col='blue',lwd=2,lty=2)

points(as.Date(pillow.dates),pillow.swe,cex=1.25,col='black',pch=16)
text(as.Date('2005-01-01'),0.95*ymax,site.name,cex=2)
legend('topright',legend=c('Pillow Obs.','ERA','PNWNAmet'),col=c('black','blue','green'),pch=16,cex=1.15)

box(which='plot')
abline(h=0)

plot(as.Date(clim.dates),clim.tas,cex=1.5,col='white',pch=16,
     xlim=c(as.Date(dst),as.Date(den)),ylim=range(c(min(clim.tasmin),max(clim.tasmax))),
     main='',xlab='Date',ylab='Temperature (degC)', cex.lab=1.25,axes=T)
axis(1,at=as.Date(c('1998-01-01','2000-01-01','2002-01-01','2004-01-01')),label=c('1998','2000','2002','2004'),cex.axis=1.25)
##axis(2,at=seq(0,yupp,round((yupp-100)/3,-2)),label=seq(0,yupp,round((yupp-100)/3,-2)),cex.axis=1.25)
lines(as.Date(clim.dates),clim.tasmax,col='red',lwd=2)
lines(as.Date(clim.dates),clim.tas,col='orange',lwd=2)
lines(as.Date(clim.dates),clim.tasmin,col='goldenrod',lwd=2)
box(which='plot')
abline(h=0)
}


plot(as.Date(pillow.dates),pillow.swe,cex=1,col='black',pch=16,
     xlim=c(as.Date('1989-08-01'),as.Date('2018-07-31')),ylim=c(0,ymax),
     main='',xlab='Date',ylab='SWE (mm)', cex.lab=1.25,axes=F)
axis(1,at=as.Date(c('1980-01-01','1990-01-01','2000-01-01','2010-01-01')),label=c('1980','1990','2000','2010'),cex.axis=1.25)
axis(2,at=seq(0,yupp,round((yupp-100)/3,-2)),label=seq(0,yupp,round((yupp-100)/3,-2)),cex.axis=1.25)

##lines(as.Date(clim.dates),swe.mean,lwd=3,col='blue')
##apply(cumulative.end.swe,2,function(y,x){lines(x,y,col='gray',lwd=1)},as.Date(clim.dates))
##lines(as.Date(clim.dates),apply(cumulative.end.swe,1,mean),col='lightblue',lwd=2)

##apply(cumulative.sim.swe,2,function(y,x){lines(x,y,col='gray',lwd=1)},as.Date(clim.dates))
##lines(as.Date(clim.dates),apply(cumulative.sim.swe,1,mean),col='darkgreen',lwd=2)
lines(as.Date(clim.dates),cumulative.sim.pr,col='green',lwd=2)
lines(as.Date(clim.dates),cumulative.end.pr,col='green',lwd=2,lty=2)
lines(as.Date(clim.dates),cumulative.snow.mean,col='blue',lwd=2)
lines(as.Date(clim.dates),cumulative.snow.end,col='blue',lwd=2,lty=2)

points(as.Date(pillow.dates),pillow.swe,cex=1,col='black',pch=16)
text(as.Date('2005-01-01'),0.95*ymax,site.name,cex=2)
legend('topright',legend=c('Pillow Obs.','ERA','PNWNAmet'),col=c('black','blue','green'),pch=16,cex=1.15)
lines(as.Date(clim.dates),results$swe*1000,col='red',lwd=2)
box(which='plot')
abline(h=0)


paper.coeffs <- list(a=-49.49,b=0.5628,c=1.5,d=1.0209)
tas <- seq(-20,20,0.05)

paper.frac <- paper.coeffs$a*(tanh(paper.coeffs$b*(tas-paper.coeffs$c))-paper.coeffs$d)
frac2 <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)

par(mar=c(5,5,2,1))
plot(tas,paper.frac,type='l',lwd=2,xlab="Daily Average Temperature (\u00B0C)",ylab='Percent Snow (%)',
              xaxs='i',yaxs='i',cex.lab=2,cex.axis=2)
lines(tas,frac2,col='red',lwd=2)
box(which='plot')


clim.subset <- as.Date(clim.dates) %in% as.Date(pillow.dates)
pillow.subset <- as.Date(pillow.dates) %in% as.Date(clim.dates)

print(mean(abs(results$swe[clim.subset]*1000 - pillow.swe[pillow.subset]),na.rm=T))
print(mean((results$swe[clim.subset]*1000 - pillow.swe[pillow.subset]),na.rm=T))


